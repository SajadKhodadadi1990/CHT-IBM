        
        forAll(phi,faceI)
        {
                    phiT[faceI] =phi[faceI];
        }
        forAll(phiT.boundaryField(),patchi)
        {
                      forAll(phiT.boundaryField()[patchi],faceI)
                      {
                      		phiT.boundaryField()[patchi][faceI]=phi.boundaryField()[patchi+1][faceI];
                      }
        }
        forAll(gammaT,cellI)
        {
                    gammaT[cellI] =gamma[cellI];
                    gammas[cellI] = 0.0;
        }
        forAll(gammaT.boundaryField(),patchi)
        { 
                      forAll(gammaT.boundaryField()[patchi],faceI)
                      {
                      		gammaT.boundaryField()[patchi][faceI]=gamma.boundaryField()[patchi+1][faceI];
                      }
        }
        DDD = gammaT * DTL + (1-gammaT) * DTS;
	fvScalarMatrix TEqn 
	(
	   fvm::ddt(T)
	 + fvm::div(phiT, T)
	 - fvm::laplacian(DDD, T)
	);
        TEqn.solve(mesh.solutionDict().solver("T"));

        forAll(Trev,cellI)
        {
                    Trev[cellI] =T[cellI];
        }
        forAll(T.boundaryField(),patchi)
        { 
                      forAll(T.boundaryField()[patchi],faceI)
                      {
                      		Trev.boundaryField()[patchi+1][faceI]=T.boundaryField()[patchi][faceI];
                      }
        }
        forAll(phiT,faceI)
        {
         	if(gamma[own[faceI]]==0&&gamma[nei[faceI]]==1)
                {

                      scalar a=DTS.value()/(mag(mesh.C()[own[faceI]]-mesh.Cf()[faceI]));
                      scalar b=DTL.value()/(mag(mesh.C()[nei[faceI]]-mesh.Cf()[faceI]));
                      scalar TS=((b/(a+b))*T[nei[faceI]]) + ((a/(a+b))*T[own[faceI]]);
                      scalar Tgrad=mag(T[nei[faceI]]-TS)/(mag(mesh.C()[nei[faceI]]-mesh.Cf()[faceI]));
                      TIbFace[nei[faceI]]=Tgrad/(TS-TRef.value());
                    //TIbFace=(T[nei[faceI]]+T[own[faceI]])/2.0;
                }
         	if(gamma[own[faceI]]==1&&gamma[nei[faceI]]==0)
                {
                      scalar a=DTS.value()/(mag(mesh.C()[nei[faceI]]-mesh.Cf()[faceI]));
                      scalar b=DTL.value()/(mag(mesh.C()[own[faceI]]-mesh.Cf()[faceI]));
                      scalar TS=((b/(a+b))*T[own[faceI]]) + ((a/(a+b))*T[nei[faceI]]);
                      scalar Tgrad=mag(T[own[faceI]]-TS)/(mag(mesh.C()[own[faceI]]-mesh.Cf()[faceI]));
                      TIbFace[own[faceI]]=Tgrad/(TS-TRef.value());
                }
        }

    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );

		nuss.open("nuss" , std::ios::out | std::ios::app);
		const labelList& ibCells = ibPatch.ibPolyPatch().ibCells();
		labelList NearCell(ibCells.size(), 0);
                 label cellll=ibCells[0];
                for(label i=1;i<ibCells.size();i++)
                {
					//Info<<"cellll:"<<newMesh.C()[ibCells[i]].component(0)<<"\t"<<newMesh.C()[cellll].component(0)<<endl;

                    if(newMesh.C()[ibCells[i]].component(0)+SMALL<newMesh.C()[cellll].component(0)||newMesh.C()[ibCells[i]].component(1)+SMALL<newMesh.C()[cellll].component(1)||newMesh.C()[ibCells[i]].component(2)+SMALL<newMesh.C()[cellll].component(2))
                    {
						cellll=ibCells[i];
					}
                }
				labelList NearCellList(ibCells.size(), 0);
				labelList checkList(ibCells.size(), 0);

				label cellI=cellll;
				label cellII=0;
				forAll(ibCells,celli)
                { 
					//~ if(celli==ibCells.size())
					//~ {
						//~ //NearCell[celli]=celli;
						//~ break;
					//~ }
					checkList[celli]=cellI;
					findCellCells(newMesh.cellCentres()[cellI],ibCells,NearCellList,newMesh.cellCentres());
						//Info<<"NearCellList["<<cellI<<"]:"<<NearCellList<<endl;
						//Info<<"checkList["<<cellI<<"]:"<<checkList<<endl;
						int JJ=1;
						bool a=true;
						bool b=true;
                        do
                        {
								bool b=true;
								//Info<<"JJ["<<cellI<<"]:"<<JJ<<endl;
								for(int i=0;i<=celli;i++)
								{
									if(NearCellList[JJ]==checkList[i])
									{
										b=false;
									}
								}
								if(b)
								{
									//Info<<"findCellCells["<<cellI<<","<<cellII<<"]:"<<NearCellList[JJ]<<endl;
									NearCell[celli]=cellI;
									cellI=NearCellList[JJ];
									a=false;
									break;
								}
							JJ++;
							if(JJ>ibCells.size())
							{
								break;
							}
						}while(a);
				}
					//Info<<"findCellCells[""]:"<<NearCell<<endl;
//--------------------------------------------------------------------------
                forAll(ibCells,celli)
                {
                  TIbFace.boundaryField()[0][celli]=TIbFace[ibCells[celli]];
                  Trev.boundaryField()[0][celli]=TIbFace[ibCells[celli]];
                }
                //----------------------------------
                 forAll(NearCell,celli)
                {
                  if(counter%200==0)
                  {
		  	nuss <<runTime.timeName()<<"\t"<<NearCell[celli]<<"\t"<<mesh.C()[NearCell[celli]].component(0)<<"\t"<<mesh.C()[NearCell[celli]].component(1)
                  	<<"\t"<<mesh.C()[NearCell[celli]].component(2)<<"\t"<<3.3*TIbFace[NearCell[celli]]<<"\n"<<endl;
                  }
                }
                //-------------------------------------
                nuss<<"\n"<<"\n"<<"\n"<<"\n"<<endl;
		nuss.close();
        }
    }
        counter++;
	Info<< "Temperature Field results = "
		   << T.weightedAverage(mesh.V()).value()
		   << "  Min(T) = " << min(T).value()
		   << "  Max(T) = " << max(T).value()
		   << endl;
