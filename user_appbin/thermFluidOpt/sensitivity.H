
Info<<"sensitivity analysis"<<endl;

fsenshMeanT=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoc;
if(DDHS)
{
 fsenshMeanT += b1*Tb*100000000.0/4180000*Q;
}
gsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);

fsenshMeanT=fsenshMeanT/gMax(mag(fsenshMeanT.primitiveField()));
gsenshPowerDiss=gsenshPowerDiss/gMax(mag(gsenshPowerDiss.primitiveField()));
gsenshVol.primitiveFieldRef()=mesh.V()/gMax(mesh.V());

#include "filter_sens.H"

gx[0]=(PowerDiss/Foam::max(PowerDissRelax-opt*0.2,PowerDissMax)-1.0);
gx[1]=V;

VecGetArray(dfdx,&dfdx_array);
VecGetArrays(dgdx,m,&dgdx_array);
for(i=0;i<n;i++)
{
 dfdx_array[i]=fsensMeanT[i]/N;//sensitivities of objective function
 dgdx_array[0][i]=gsensPowerDiss[i]/N;//sensitivities of constraint function
 dgdx_array[1][i]=gsensVol[i]/N;
}
VecRestoreArray(dfdx,&dfdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

Info<< "\n run MMA \n" << endl;
mma->SetOuterMovelimit(Xmin,Xmax,movlim,xmma,xmin,xmax);
mma->Update(xmma,dfdx,gx,dgdx,xmin,xmax,raa0,mma_dec,mma_inc);

VecGetArray(xmma,&xmma_array);
for(i=0;i<n;i++)
{
   x[i]=xmma_array[i];
}
VecRestoreArray(xmma,&xmma_array);

#include "filter_x.H"
