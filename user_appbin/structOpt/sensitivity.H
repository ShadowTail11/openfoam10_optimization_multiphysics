
Info<<"sensitivity analysis"<<endl;
fsensh=-gradD&&(rho*(3*xh*xh*(E-Emin)/(2.0*(1.0 + Po))*twoSymm(gradD) + (3*xh*xh*Po*(E-Emin)/((1.0 + Po)*(1.0 - 2.0*Po))*I)*tr(gradD)));  
for(i=0;i<n;i++)
{
 fsensh[i]=fsensh[i]*mesh.V()[i];
}
fsensh=fsensh/gMax(mag(fsensh.primitiveField()));
gsensh.primitiveFieldRef()=mesh.V()/gMax(mesh.V());

#include "filter_chainrule.H"

gx[0]=V;//constraint function

VecGetArray(dfdx,&dfdx_array);
VecGetArrays(dgdx,m,&dgdx_array);
for(i=0;i<n;i++)
{
 dfdx_array[i]=fsens[i]/N;//sensitivities of objective function
 dgdx_array[0][i]=gsens[i]/N;//sensitivities of constraint function
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
