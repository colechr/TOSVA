//! \file  clock.c
//! \brief Contains simulation time function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

PetscErrorCode adjustTimeStep (domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;
    flags_   *flags   = domain[0].access.flags;

    PetscInt  flag       = 0;
    PetscReal cfl        = 1e10;
    PetscReal maxU       = 0.0;
    PetscReal dxByU_min  = 1e10;
    PetscReal dxBy__min  = 1e10;
    PetscReal dxBygTau_min  = 1e10;
    PetscReal dxByaTau_min  = 1e10;
    cellIds   maxUCell;

    clock_        *clock = domain[0].clock;

    // save old time step
    clock->dtOld = clock->dt;

    for(PetscInt d=0; d<nDomains; d++)
    {
        acquisition_  *acquisition = domain->acquisition;

        if(flags->isScalarMomentsActive)
        {
            // time step imposed by the flow on this domain
            timeStepInfo(&domain[d], clock, dxByU_min, maxU, maxUCell);

            if (flags->isSediFluxActive && (clock->it != clock->itStart)) //skip first iteration for stability, Allows U-field to develop first.
            {
                timeStepInfoSMSed(&domain[d], clock, dxBygTau_min);

                if (dxBygTau_min == 1e10)
                {
                    char error[512];
                    sprintf(error, "DIVERGED SMs. EIther Tau is 0 everywhere, or nan exists. Try lower CFL or lower GSD/GMD");
                    fatalErrorInFunction("adjustTimeStep", error);
                }
            }

            if (flags->isDeviFluxActive && (clock->it != clock->itStart)) //skip first iteration for stability; acceleration term needed but not yet availabile until later iterations.
            {
                timeStepInfoSMDev(&domain[d], clock, dxByaTau_min);

                if (dxByaTau_min == 1e10)
                {
                    char error[512];
                    sprintf(error, "DIVERGED SMs. EIther Tau is 0 everywhere, or nan exists. Try lower CFL or lower GSD/GMD");
                    fatalErrorInFunction("adjustTimeStep", error);
                }
            }
        }
        else
        {
            // set cfl
            cfl = PetscMin(cfl, clock->cfl);

            // time step imposed by the flow on this domain
            timeStepInfo(&domain[d], clock, dxByU_min, maxU, maxUCell);
        }

        // time step imposed by the precursor on this domain (must be corrected to syncronize time steps)
        if(domain[d].flags.isConcurrentPrecursorActive)
        {
            //timeStepInfo(domain[d].abl->precursor->domain, clock, dxByU_min, maxU, maxUCell);
        }

        // try to guess a uniform predicted time step for the acquisition
        PetscReal predictedDt;

        // output fields
        if(flags->isAdjustableTime)
        {
            //take minimum from covective (divergence), sedimentive, and deviative CFL. Takes convective if SM is off.
            dxBy__min = PetscMin(dxByU_min, PetscMin(dxBygTau_min, dxByaTau_min));

            if(flags->isScalarMomentsActive)
            {
                if (dxBy__min == dxByU_min)
                {
                    clock->dt   = clock->cfl * dxBy__min;
                    cfl = PetscMin(cfl, clock->cfl);
                }
                else
                {
                    clock->dt   = clock->cflSM * dxBy__min;
                    cfl = PetscMin(cfl, clock->cflSM);
                }
            }
            else
            {
                // 2. takes the local ratio
                clock->dt   = clock->cfl * dxBy__min;
            }


            PetscReal timeStart;
            PetscReal timeInterval;

            timeStart    = clock->startTime;
            timeInterval = domain[d].io->timeInterval;

            timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
            predictedDt  = currentDistanceToWriteTime(clock, timeStart, timeInterval);

            // averaged tke
            if(domain[d].io->TKE)
            {
                timeStart    = domain[d].io->tkeStartTime;
                timeInterval = domain[d].io->tkePrd;

                timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            // averaged fields
            if(domain[d].io->averaging)
            {
                timeStart    = domain[d].io->avgStartTime;
                timeInterval = domain[d].io->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            // phase averaged fields
            if(domain[d].io->phaseAveraging)
            {
                timeStart    = domain[d].io->phAvgStartTime;
                timeInterval = domain[d].io->phAvgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            // ke budgets
            if(domain[d].io->keBudgets)
            {
                timeStart    = acquisition->keBudFields->avgStartTime;
                timeInterval = acquisition->keBudFields->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            if(flags->isIBMActive)
            {
                ibm_ *ibm = domain[d].ibm;

                if(domain[d].ibm->intervalType == "adjustableTime")
                {
                    timeStart    = ibm->timeStart;
                    timeInterval = ibm->timeInterval;

                    timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }
            }

            if(flags->isWindFarmActive)
            {
                if(domain[d].farm->intervalType == "adjustableTime")
                {
                    timeStart    = domain[d].farm->timeStart;
                    timeInterval = domain[d].farm->timeInterval;

                    timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }
            }

            // catalyst
            if(flags->isPvCatalystActive)
            {
                #ifdef USE_CATALYST
                if(domain[d].io->outputTypeCatalyst == "adjustableTime")
                {
                    timeStart    = domain[d].io->startTimeCatalyst;
                    timeInterval = domain[d].io->timeIntervalCatalyst;

                    timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }
                #endif
            }

            if(flags->isAquisitionActive)
            {
                // 3LM averaging for background domain only
                if(acquisition->isAverage3LMActive)
                {
                    timeStart    = domain[0].acquisition->LM3->avgStartTime;
                    timeInterval = domain[0].acquisition->LM3->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }

                // ABL perturbations
                if(acquisition->isPerturbABLActive)
                {
                    timeStart    = domain[0].acquisition->perturbABL->avgStartTime;
                    timeInterval = domain[0].acquisition->perturbABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }

                // ABL averaging for background domain only
                if(acquisition->isAverageABLActive)
                {
                    timeStart    = domain[0].acquisition->statisticsABL->avgStartTime;
                    timeInterval = domain[0].acquisition->statisticsABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }

                if(acquisition->isProbesActive)
                {
                    for(PetscInt r=0; r<acquisition->probes->nRakes; r++)
                    {
                        if(acquisition->probes->rakes[r].intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->probes->rakes[r].timeStart;
                            timeInterval = acquisition->probes->rakes[r].timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));

                            if(acquisition->probes->allSameIO) break;
                        }
                    }
                }

                if(acquisition->isSectionsActive)
                {
                    // iSections
                    if(acquisition->iSections->available)
                    {
                        if(acquisition->iSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->iSections->timeStart;
                            timeInterval = acquisition->iSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                        }
                    }

                    // iSections
                    if(acquisition->jSections->available)
                    {
                        if(acquisition->jSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->jSections->timeStart;
                            timeInterval = acquisition->jSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                        }
                    }

                    // kSections
                    if(acquisition->kSections->available)
                    {
                        if(acquisition->kSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->kSections->timeStart;
                            timeInterval = acquisition->kSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                        }
                    }
                }
            }

            // concurrent precursor
            if(domain[d].flags.isConcurrentPrecursorActive)
            {
                /*
                if(domain[d].abl->precursor->domain->flags.isAquisitionActive)
                {
                    if(domain[d].abl->precursor->domain->acquisition->isAverageABLActive)
                    {
                        timeStart    = domain[d].abl->precursor->domain->acquisition->statisticsABL->avgStartTime;
                        timeInterval = domain[d].abl->precursor->domain->acquisition->statisticsABL->avgPrd;

                        timeStepSet(clock, timeStart, timeInterval, dxBy__min, flag, cfl);
                        predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                    }
                }
                */
            }

            // set time step as the gcd of all constraints
            if(clock->it == 0)
            {
                // scale max uniform dt due to acquisition so that it complies the CFL
                while(predictedDt / dxBy__min > clock->cfl)
                {
                    predictedDt /= 2.0;
                }

                // save
                clock->acquisitionDt = predictedDt;
            }

            // check if must limit dt due to ALM rotation
            if(flags->isWindFarmActive)
            {
                farm_ *farm = domain[d].farm;

                // set checkCFL flag
                if(clock->it == 0)
                {
                    // loop over each wind turbine
                    for(PetscInt t=0; t<farm->size; t++)
                    {
                        if((*farm->turbineModels[t]) == "ALM")
                        {
                            farm->checkCFL = 1;
                        }
                    }
                }

                if(farm->checkCFL)
                {
                    computeMaxTipSpeed(domain[d].farm);

                    PetscReal dtFarm = clock->dxMin / farm->maxTipSpeed;

                    clock->dt    = std::min(clock->dt, dtFarm);
                }
            }

            // check cfl limit for rotating ibm
            if(flags->isIBMActive)
            {
                ibm_ *ibm = domain[d].ibm;

                for (PetscInt i=0; i < ibm->numBodies; i++)
                {
                    ibmObject   *ibmBody = ibm->ibmBody[i];

                    if(ibm->dynamic)
                    {
                        if(ibmBody->bodyMotion == "rotation")
                        {
                            ibmRotation *ibmRot  = ibmBody->ibmRot;

                            PetscReal maxSpeed   = ibmRot->maxR * ibmRot->angSpeed;

                            PetscReal dtIBM      = dxBy__min*maxU / maxSpeed;

                            clock->dt    = std::min(clock->dt, dtIBM);
                        }
                    }
                }
            }
        }
        else
        {
            if(clock->it>clock->itStart) clock->cfl =  clock->dt / dxBy__min;

            // added by Arjun for fixed time step as a last control
            /*
			if(clock->cfl > 1.0)
            {
                clock->dt = clock->dt/2.0;
            }
			*/
        }
    }

    // prevent time step from increasing too fast
    if(flags->isAdjustableTime)
    {
        if(clock->dt > 1.5 * clock->dtOld)
        {
            clock->dt    = 1.5 * clock->dtOld;
        }
    }

    // avoid too small time step (will not write, have to fix this)
    if (clock->dt < 1e-10)
    {
        clock->dt = clock->startDt;
    }

    // discard all previous changes if time step is fixed as acquistion gcd
    if(flags->isAdjustableTime == 2)
    {
        if(clock->acquisitionDt / dxBy__min > 1.0)
        {
            clock->acquisitionDt /= 2.0;
        }

        clock->dt = clock->acquisitionDt;
    }

    // make sure to hit last time
    if(clock->time + clock->dt > clock->endTime)
    {
        clock->dt = clock->endTime - clock->time;
    }

    if(flags->isAdjustableTime)
    {
        clock->time = clock->time + clock->dt;
    }
    else
    {
        //to exit when reaching end time
        if(clock->time + clock->dt > clock->endTime)
        {
            clock->time = clock->time + clock->dt;
        }
        else
        {
            // this ensures there is no floating point addition error
            clock->time = clock->startTime + (clock->it + 1) * clock->dt;
        }
    }

    if(!flags->isScalarMomentsActive)
    {
        cfl = clock->dt / dxBy__min;
    }

    const char *cfl_label = "DIV";

    if (dxBy__min == dxBygTau_min)
    {
        cfl_label = "SED";
        PetscPrintf(PETSC_COMM_WORLD, "\n\nWARNING ... Scalar moments sedFlux is limiting momentum time step! Lower GMD or GSD, or turn sedFlux off to avoid.\n\n");
    }
    else if (dxBy__min == dxByaTau_min)
    {
        cfl_label = "DEV";
        PetscPrintf(PETSC_COMM_WORLD, "\n\nWARNING ... Scalar moments devFlux is limiting momentum time step! Lower GMD or GSD, or turn devFlux off to avoid.\n\n");
    }


    PetscPrintf(PETSC_COMM_WORLD, "\n\nTime: %lf\n\n", clock->time);

    if(clock->it==clock->itStart)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Iteration = %ld, CFL_%s = %lf, uMax = %.6f, dt = %.6f, dtMaxCFL = %.6f, adjust due to write flag: %ld\n", clock->it, cfl_label, cfl, maxU, clock->dt, dxBy__min, flag);
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "Iteration = %ld, CFL_%s = %lf, uMax = %.6f (i,j,k = %ld, %ld, %ld), dt = %.6f, dtMaxCFL = %.6f, adjust due to write flag: %ld\n", clock->it, cfl_label, cfl, maxU, maxUCell.i, maxUCell.j, maxUCell.k, clock->dt, dxBy__min, flag);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode timeStepInfo(domain_ *domain, clock_ *clock, PetscReal &dxByU_min, PetscReal &maxU, cellIds &maxUCell)
{
    flags_        *flags = domain->access.flags;
    mesh_         *mesh  = domain->mesh;
    ueqn_         *ueqn  = domain->ueqn;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucat, ***ucont;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert, ***aj;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    PetscReal lmaxU      = 0.0,  gmaxU      = 0.0;
    PetscReal ldx        = 1e10;
    PetscReal ldi_min    = 1e10, ldj_min    = 1e10, ldk_min = 1e10;
    PetscReal ldxByU_min = 1e10, gdxByU_min = 1e10;
    cellIds   lmaxUCell; lmaxUCell.i = 0; lmaxUCell.j = 0; lmaxUCell.k = 0;
    cellIds   gmaxUCell; gmaxUCell.i = 0; gmaxUCell.j = 0; gmaxUCell.k = 0;

    DMDAVecGetArray(fda, ueqn->Ucat,   &ucat);
    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    // time step due to flow restrictions
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // cartesian velocity magnitude
                PetscReal Umag = nMag(ucat[k][j][i]);

                // maximum overall velocity in this proc
                if(Umag > lmaxU)
                {
                    lmaxU = Umag;
                    lmaxUCell.i = i;
                    lmaxUCell.j = j;
                    lmaxUCell.k = k;
                }

                // compute cell sizes
                PetscReal ldi = 1./aj[k][j][i]/nMag(csi[k][j][i]);
                PetscReal ldj = 1./aj[k][j][i]/nMag(eta[k][j][i]);
                PetscReal ldk = 1./aj[k][j][i]/nMag(zet[k][j][i]);

                // minimum overall GCC sizes in this proc
                if(clock->it == clock->itStart)
                {
                    ldi_min = PetscMin(ldi_min, ldi);
                    ldj_min = PetscMin(ldj_min, ldj);
                    ldk_min = PetscMin(ldk_min, ldk);
                }

                // compute directional velocity in curvilinear coordinates
                PetscReal Vi = fabs(0.5 * (ucont[k][j][i].x + ucont[k][j][i-1].x) / nMag(csi[k][j][i])),
                          Vj = fabs(0.5 * (ucont[k][j][i].y + ucont[k][j-1][i].y) / nMag(eta[k][j][i])),
                          Vk = fabs(0.5 * (ucont[k][j][i].z + ucont[k-1][j][i].z) / nMag(zet[k][j][i]));

                // compute dxByU
                PetscReal ldxByU = PetscMin(ldi/Vi, PetscMin(ldj/Vj, ldk/Vk));

                if(isFluidCell(k, j, i, nvert))
                {
                    // min overall cell size
                    if(clock->it == clock->itStart)
                    {
                        ldx = PetscMin(ldi_min, PetscMin(ldj_min, ldk_min));
                    }

                    // min overall dxByU
                    ldxByU_min = PetscMin(ldxByU_min, ldxByU);
                }
            }
        }
    }

    if(clock->it == clock->itStart)
    {
        ldx = PetscMin(ldx, clock->dxMin);
        MPI_Allreduce(&ldx, &clock->dxMin, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);
    }

    // max in this domain
    MPI_Allreduce(&lmaxU,        &gmaxU,      1, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);
    MPI_Allreduce(&ldxByU_min,   &gdxByU_min, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

    // check if maxU is in this processor and print where it is
    if(lmaxU != gmaxU)
    {
        lmaxUCell.i = 0; lmaxUCell.j = 0; lmaxUCell.k = 0;
    }

    MPI_Allreduce(&lmaxUCell, &gmaxUCell, 3, MPIU_INT, MPIU_MAX, mesh->MESH_COMM);

    // max in all domains and save ids
    if(gmaxU > maxU)
    {
        maxU     = gmaxU;
        maxUCell = gmaxUCell;
    }
    if(gdxByU_min < dxByU_min)
    {
        dxByU_min = gdxByU_min;
    }

    DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode timeStepInfoSMSed(domain_ *domain, clock_ *clock, PetscReal &dxBygTau_min)
{
    flags_        *flags = domain->access.flags;
    mesh_         *mesh  = domain->mesh;
    ueqn_         *ueqn  = domain->ueqn;
    SMObj_        *smObject  = domain->smObject;

    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucat, ***ucont;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert, ***aj;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***TauP0, ***TauP1, ***TauP2;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    PetscReal ldx        = 1e10;
    PetscReal ldi_min    = 1e10, ldj_min    = 1e10, ldk_min = 1e10;
    PetscReal ldxBygTau_min = 1e10, gdxBygTau_min = 1e10;

    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    // time step due to particle sedimentiation
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(!isIBMSolidCell(k, j, i, nvert))
                {

                    // compute cell sizes
                    PetscReal ldi = 1./aj[k][j][i]/nMag(csi[k][j][i]);
                    PetscReal ldj = 1./aj[k][j][i]/nMag(eta[k][j][i]);
                    PetscReal ldk = 1./aj[k][j][i]/nMag(zet[k][j][i]);

                    // compute dxBygTau
                    PetscReal ldxBygTau0 = PetscMin(ldi/(9.81*TauP0[k][j][i]), PetscMin(ldj/(9.81*TauP0[k][j][i]), ldk/(9.81*TauP0[k][j][i])));
                    PetscReal ldxBygTau1 = PetscMin(ldi/(9.81*TauP1[k][j][i]), PetscMin(ldj/(9.81*TauP1[k][j][i]), ldk/(9.81*TauP1[k][j][i])));
                    PetscReal ldxBygTau2 = PetscMin(ldi/(9.81*TauP2[k][j][i]), PetscMin(ldj/(9.81*TauP2[k][j][i]), ldk/(9.81*TauP2[k][j][i])));

                    PetscReal ldxBygTau = PetscMin(ldxBygTau0, PetscMin(ldxBygTau1, ldxBygTau2));

                    // min overall dxByU
                    ldxBygTau_min = PetscMin(ldxBygTau_min, ldxBygTau);

                }
            }
        }
    }

    // max in this domain
    MPI_Allreduce(&ldxBygTau_min,   &gdxBygTau_min, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

    if(gdxBygTau_min < dxBygTau_min)
    {
        dxBygTau_min = gdxBygTau_min;
    }

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode timeStepInfoSMDev(domain_ *domain, clock_ *clock, PetscReal &dxByaTau_min)
{
    flags_        *flags = domain->access.flags;
    mesh_         *mesh  = domain->mesh;
    ueqn_         *ueqn  = domain->ueqn;
    SMObj_        *smObject  = domain->smObject;

    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert, ***aj;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***TauP0, ***TauP1, ***TauP2;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    PetscReal     du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

    PetscReal     dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz;
    PetscReal     dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords

    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc;      // surface area vectors components

    Cmpnts	      ***icsi, ***ieta, ***izet;
    Cmpnts	      ***jcsi, ***jeta, ***jzet;
    Cmpnts	      ***kcsi, ***keta, ***kzet;

    PetscReal     ***iaj, ***jaj, ***kaj;

    PetscInt      ***markVent;

    Cmpnts        ***ucat, ***ucat_o, matDerI, matDerJ, matDerK, matDer;

    PetscReal ldx        = 1e10;
    PetscReal ldxByaTau_min = 1e10, gdxByaTau_min = 1e10;

    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, ueqn->Ucat_o, &ucat_o);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);

    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);

    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
    DMDAVecGetArray(da, mesh->lJAj, &jaj);
    DMDAVecGetArray(da, mesh->lKAj, &kaj);

    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    // time step due to particle deviation restrictions
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(!isIBMSolidCell(k, j, i, nvert))
                {
                    //i-faces
                    csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                    eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                    zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
                    ajc  = iaj[k][j][i];

                    // compute cartesian velocity derivatives w.r.t. curvilinear coords
                    Compute_du_i
                    (   mesh, i, j, k, mx, my, mz, ucat, nvert,
                        &dudc, &dvdc, &dwdc,
                        &dude, &dvde, &dwde,
                        &dudz, &dvdz, &dwdz
                    );

                    // compute cartesian velocity derivatives w.r.t cartesian coords
                    Compute_du_dxyz
                    (
                        mesh,
                        csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                        zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                        dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                        &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                    );

                    matDerI.x = (ucat[k][j][i].x - ucat_o[k][j][i].x)/clock->dt + ucat[k][j][i].x * du_dx + ucat[k][j][i].y * du_dy + ucat[k][j][i].z * du_dz;
                    matDerI.y = (ucat[k][j][i].y - ucat_o[k][j][i].y)/clock->dt + ucat[k][j][i].x * dv_dx + ucat[k][j][i].y * dv_dy + ucat[k][j][i].z * dv_dz;
                    matDerI.z = (ucat[k][j][i].z - ucat_o[k][j][i].z)/clock->dt + ucat[k][j][i].x * dw_dx + ucat[k][j][i].y * dw_dy + ucat[k][j][i].z * dw_dz;

                    //j-face
                    csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                    eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                    zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
                    ajc  = jaj[k][j][i];

                    // compute cartesian velocity derivatives w.r.t. curvilinear coords
                    Compute_du_j
                    (   mesh, i, j, k, mx, my, mz, ucat, nvert,
                        &dudc, &dvdc, &dwdc,
                        &dude, &dvde, &dwde,
                        &dudz, &dvdz, &dwdz
                    );

                    // compute cartesian velocity derivatives w.r.t cartesian coords
                    Compute_du_dxyz
                    (
                        mesh,
                        csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                        zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                        dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                        &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                    );

                    matDerJ.x = (ucat[k][j][i].x - ucat_o[k][j][i].x)/clock->dt + ucat[k][j][i].x * du_dx + ucat[k][j][i].y * du_dy + ucat[k][j][i].z * du_dz;
                    matDerJ.y = (ucat[k][j][i].y - ucat_o[k][j][i].y)/clock->dt + ucat[k][j][i].x * dv_dx + ucat[k][j][i].y * dv_dy + ucat[k][j][i].z * dv_dz;
                    matDerJ.z = (ucat[k][j][i].z - ucat_o[k][j][i].z)/clock->dt + ucat[k][j][i].x * dw_dx + ucat[k][j][i].y * dw_dy + ucat[k][j][i].z * dw_dz;

                    //k-face
                    csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                    eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                    zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
                    ajc  = kaj[k][j][i];

                    // compute cartesian velocity derivatives w.r.t. curvilinear coords
                    Compute_du_k
                    (   mesh, i, j, k, mx, my, mz, ucat, nvert,
                        &dudc, &dvdc, &dwdc,
                        &dude, &dvde, &dwde,
                        &dudz, &dvdz, &dwdz
                    );

                    // compute cartesian velocity derivatives w.r.t cartesian coords
                    Compute_du_dxyz
                    (
                        mesh,
                        csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                        zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                        dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                        &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                    );

                    matDerK.x = (ucat[k][j][i].x - ucat_o[k][j][i].x)/clock->dt + ucat[k][j][i].x * du_dx + ucat[k][j][i].y * du_dy + ucat[k][j][i].z * du_dz;
                    matDerK.y = (ucat[k][j][i].y - ucat_o[k][j][i].y)/clock->dt + ucat[k][j][i].x * dv_dx + ucat[k][j][i].y * dv_dy + ucat[k][j][i].z * dv_dz;
                    matDerK.z = (ucat[k][j][i].z - ucat_o[k][j][i].z)/clock->dt + ucat[k][j][i].x * dw_dx + ucat[k][j][i].y * dw_dy + ucat[k][j][i].z * dw_dz;

                    //printf("%f %f %f %f %f %f\n", ucat[k][j][i].z, ucat_o[k][j][i].z, clock->dt, dw_dx, dw_dy, dw_dz);

                    matDer.x = PetscMax(fabs(matDerI.x), PetscMax(fabs(matDerK.x), fabs(matDerJ.x)));
                    matDer.y = PetscMax(fabs(matDerI.y), PetscMax(fabs(matDerK.y), fabs(matDerJ.y)));
                    matDer.z = PetscMax(fabs(matDerI.z), PetscMax(fabs(matDerK.z), fabs(matDerJ.z)));

                    // compute cell sizes
                    PetscReal ldi = 1./aj[k][j][i]/nMag(csi[k][j][i]);
                    PetscReal ldj = 1./aj[k][j][i]/nMag(eta[k][j][i]);
                    PetscReal ldk = 1./aj[k][j][i]/nMag(zet[k][j][i]);

                    // compute dxBygTau
                    PetscReal ldxByaTau0 = PetscMin(ldi/(matDer.x*TauP0[k][j][i]), PetscMin(ldj/(matDer.y*TauP0[k][j][i]), ldk/(matDer.z*TauP0[k][j][i])));
                    PetscReal ldxByaTau1 = PetscMin(ldi/(matDer.x*TauP1[k][j][i]), PetscMin(ldj/(matDer.y*TauP1[k][j][i]), ldk/(matDer.z*TauP1[k][j][i])));
                    PetscReal ldxByaTau2 = PetscMin(ldi/(matDer.x*TauP2[k][j][i]), PetscMin(ldj/(matDer.y*TauP2[k][j][i]), ldk/(matDer.z*TauP2[k][j][i])));

                    PetscReal ldxByaTau = PetscMin(ldxByaTau0, PetscMin(ldxByaTau1, ldxByaTau2));

                    // min overall dxByU
                    ldxByaTau_min = PetscMin(ldxByaTau_min, ldxByaTau);

                }
            }
        }
    }

    // max in this domain
    MPI_Allreduce(&ldxByaTau_min,   &gdxByaTau_min, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

    if(gdxByaTau_min < dxByaTau_min)
    {
        dxByaTau_min = gdxByaTau_min;
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat_o, &ucat_o);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    return(0);
}
