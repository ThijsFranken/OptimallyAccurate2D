%_______SINC INTERPOLATION FOR 3PT, 5PT AND 7PT HETEROGENEOUS SCHEME_______
%__________________________________________________________________________


%----defining parameters----
xm = 0; 
zm = 0;

  
dx = 1.d0;
dz = 1.d0;
ndis = 100;
smalldx = dx/ndis;
smalldz = dz/ndis;
dx2 = dx * dx;
dz2 = dz * dz;
dxdz = dx * dz;



mx=2;
mz=2;

nx=0;
nz=0;

pi=3.141592653589793238;
%xm1 = xm + deltax;
                             % for discretizing steps

trialfunction=1; % 1 for linear spline, 2 for sinc
npTF=3;          % number of points in scheme (3, 5, 7)
ngrid=(npTF-1)/2;   %Boundary rows

mxmin = 1;
mxmax = npTF;
mzmin = 1;
mzmax = npTF; 
% IX = 5; % the number of nodes in x
% IZ = 5; %tge number of nodes in z

% definition of T and H matrix

T=0;
% H1=zeros(mx,mz,nx,nz);
% H2=zeros(mx,mz,nx,nz);

% n=0;
% m=0;
% MN1=zeros(n,m); %Matrix describing integration range along boundaries
% d=0;            %Kronecker Delta: required to solve a numerical inconsistency
                %between diaganoal and non-diagonal matrix elements, because of
                %a range variation between these cases.


% tlength= 100;   %window length
% iomega=100;
% delta_omega=(2*pi)/tlength;
% iomegamin=0;
% iomegamax=iomega*delta_omega;
% 
% k=0;            %wavenumber = 0 for 1D case
% g=zeros(NX,1);
% B=zeros(NX,NX);
% 
% Source_Pos = (NX+1)/2; % Source position [m]
% mindist = NX;

% for i = 2:NX
%     
%     dist = Source_Pos - deltax*i ;
%     
%     if abs(dist) < mindist
%     mindist = dist ;
%     isource = i ;
%     end
%     
% end
    


%_______________________TRIAL FUNCTIONS_______________________

%Spline interpolation:
if trialfunction == 1;                      %enter 1 for spline interpolation
    for i = -ngrid*ndis:0;                            %for xm-1 < x < xm
        
        x = xm + (i/ndis)*dx;
        phix(i+1+ndis)=(x+dx)/dx;
        phixderiv(i+1+ndis)=1./dx;
        
        z = zm + (i/ndis)*dz;
        phiz(i+1+ndis)=(z+dz)/dz;
        phizderiv(i+1+ndis)=1./dz;
    end

    for i = 0:ngrid*ndis;                             %for xm < x < xm+1
        
        x = xm + (i/ndis)*dx;
        phix(i+1+ndis)=(-x+dx)/dx;
        phixderiv(i+1+ndis)=-1./dx;
        
        z = zm + (i/ndis)*dz;
        phiz(i+1+ndis)=(-z+dz)/dz;
        phizderiv(i+1+ndis)=-1./dz;
    end

%Sinc interpolation:
elseif trialfunction == 2;                  %enter 2 for sinc interpolation
    
        
        for i = -ngrid*ndis:ngrid*ndis;           %range depends on npTF; points in scheme
            x = xm+(i/ndis)*dx;
            xx= (x-xm)/dx;
            phix(i+1+ngrid*ndis) = sinc(xx);
            
            z = zm+(i/ndis)*dz;
            zz= (z-zm)/dz;
            phiz(i+1+ngrid*ndis) = sinc(zz);
            
            if i == 0;
                   phixderiv(i+1+ngrid*ndis)=0;
                   phizderiv(i+1+ngrid*ndis)=0;
            else
            %phimderiv(i+1+ndis) = (phim(i+2+ndis)-phim(i+1+ndis))/(deltax/ndis);
                    phixderiv(i+1+ngrid*ndis) = pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx); %Analytical solution
                    phizderiv(i+1+ngrid*ndis) = pi*cos(pi*zz)/(pi*zz) - sin(pi*zz)/(pi*zz*zz);
            end
        end
end


%_______________________ELASTIC PROPERTIES________________________


% linearly increasing rho: rho=a_rho*X +b_rho
a_rho=0.;         % a_rho=0 for homogeneous
b_rho=1.;

% linearly increasing mu: mu=a_mu*X +b_mu
a_mu=0.;          % a_mu=0 for homogeneous
b_mu=1.;

%____________________COMPUTATION OF T AND H MATRICES______________________
% 


mx=2;
mz=2;

markers = zeros(2*ngrid*ndis+1,2*ngrid*ndis+1,mxmax-mxmin+1,mzmax-mzmin+1);
% mx, mz centered on 2. Integrated over nx, nz:

for nx = mxmin:mxmax;
    for nz = mxmin:mxmax;
        
        T(mx,mz,nx,nz)=0;
        H11(mx,mz,nx,nz)=0;
        H13(mx,mz,nx,nz)=0;
        H31(mx,mz,nx,nz)=0;
        H33(mx,mz,nx,nz)=0;
        
        for inx = -ngrid*ndis:ngrid*ndis;
            imx = inx-(mx-nx)*ndis;
            
              if ((imx-(-ngrid*ndis))*(imx-(ngrid*ndis-1)))<=0
                
                   for inz = -ngrid*ndis:ngrid*ndis;
                        imz = inz-(mz-nz)*ndis;
                        
                        if ((imz-(-ngrid*ndis))*(imz-(ngrid*ndis-1)))<=0
                            markers(imx+ngrid*ndis+1,imz+ngrid*ndis+1,nx,nz)=1; 
                            T(mx,mz,nx,nz) = T(mx,mz,nx,nz) + phix(imx+1+ngrid*ndis)*phix(inx+1+ngrid*ndis)*phiz(imz+1+ngrid*ndis)*phiz(inz+1+ngrid*ndis)*smalldx*smalldz;
                            H11(mx,mz,nx,nz) = H11(mx,mz,nx,nz) + phixderiv(imx+1+ngrid*ndis)*phixderiv(inx+1+ngrid*ndis)*phiz(imz+1+ngrid*ndis)*phiz(inz+1+ngrid*ndis)*smalldx*smalldz;
                            H13(mx,mz,nx,nz) = H13(mx,mz,nx,nz) + phixderiv(imx+1+ngrid*ndis)*phix(inx+1+ngrid*ndis)*phiz(imz+1+ngrid*ndis)*phizderiv(inz+1+ngrid*ndis)*smalldx*smalldz;
                            H31(mx,mz,nx,nz) = H31(mx,mz,nx,nz) + phix(imx+1+ngrid*ndis)*phixderiv(inx+1+ngrid*ndis)*phizderiv(imz+1+ngrid*ndis)*phiz(inz+1+ngrid*ndis)*smalldx*smalldz;
                            H33(mx,mz,nx,nz) = H33(mx,mz,nx,nz) + phix(imx+1+ngrid*ndis)*phix(inx+1+ngrid*ndis)*phizderiv(imz+1+ngrid*ndis)*phizderiv(inz+1+ngrid*ndis)*smalldx*smalldz;
                        end
                   end
              end
        end
    end
end

pcolor(markers(:,:,1,1));
% T11 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(1*ndis)+phix(1*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(1*ndis)+phiz(1*ndis+1))*smalldx*smalldz
% T12 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(1*ndis)+phix(1*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*smalldx*smalldz
% T13 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(1*ndis)+phix(1*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(3*ndis)+phiz(3*ndis+1))*smalldx*smalldz
% % 
% T21 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(2*ndis)+phix(2*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(1*ndis)+phiz(1*ndis+1))*smalldx*smalldz
% T22 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(2*ndis)+phix(2*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*smalldx*smalldz
% T23 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(2*ndis)+phix(2*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(3*ndis)+phiz(3*ndis+1))*smalldx*smalldz
% % 
% T31 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(3*ndis)+phix(3*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(1*ndis)+phiz(1*ndis+1))*smalldx*smalldz
% T32 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(3*ndis)+phix(3*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*smalldx*smalldz
% T33 = T+ (phix(2*ndis)+phix(2*ndis+1))*(phix(3*ndis)+phix(3*ndis+1))*(phiz(2*ndis)+phiz(2*ndis+1))*(phiz(3*ndis)+phiz(3*ndis+1))*smalldx*smalldz




% %Computation of bulk matrix, excluding BR:
% for MX = 2
%     for MZ = 2
%         for NX = 1:npTF
%             for NZ = 1:npTF 
%     
% %     xm=deltax*(m-1);
% %     x=xm+(i/ndis)*deltax;
% %     rho=a_rho*x+b_rho;
% %     mu=a_mu*x+b_mu;
%     
%                 for r=-ngrid:ngrid;               %range of nodes from central node
%         
%         for n=m+r;              %computation of row elements
%         
%             T(m,n)=0.;
%             H1(m,n)=0.;
%             H2(m,n)=0.;
%             H(m,n)=0.;
%             
%             for i=-(ngrid-abs(r))*ndis:ngrid*ndis-1;      %integration range
%             
%                 T(m,n)=T(m,n) + phix(i+1+(ngrid-abs(r))*ndis)*rho*phix(i+1+ngrid*ndis)*deltax/ndis;
%                 H1(m,n)=H1(m,n) + phix(i+1+(ngrid-abs(r))*ndis)*mu*phix(i+1+ngrid*ndis)*deltax/ndis;
%                 H2(m,n)=H2(m,n) + phixderiv(i+1+(ngrid-abs(r))*ndis)*mu*phixderiv(i+1+ngrid*ndis)*deltax/ndis;
%                 H(m,n)=(k^2)*H1(m,n)+H2(m,n);
%             end
%         end
%     end
% end

%Computation of BR nodes:

% %First BR nodes
% for m=1:ngrid;
%     
%     xm=deltax*(m-1);
%     x=xm+(i/ndis)*deltax;
%     rho=a_rho*x+b_rho;
%     mu=a_mu*x+b_mu;
%     
%     for r=1-m:ngrid;               %range of nodes corrected for boundary
%         
%         for n=m+r;              %computation of row elements
%             
%             T(m,n)=0.;
%             H1(m,n)=0.;
%             H2(m,n)=0.;
%             H(m,n)=0.;
%             MN1 = [0 0 0;0 -1 -1;0 -1 -2];            %matrix describing integration range along lower boundary
%             
%                  if m==n;                             %Kronecker Delta
%                     d=1;
%                  elseif m~=n;
%                         d=0;
%                  end
%             
%             if n >= ngrid+1;                               %Known element values within boundary row;
%                                                         %follows the conventional computation for bulk matrix
%                 for i=-(ngrid-abs(r))*ndis:ngrid*ndis-1;      %integration range
%             
%                 T(m,n)=T(m,n) + phix(i+1+(ngrid-abs(r))*ndis)*rho*phix(i+1+ngrid*ndis)*deltax/ndis;
%                 H1(m,n)=H1(m,n) + phix(i+1+(ngrid-abs(r))*ndis)*mu*phix(i+1+ngrid*ndis)*deltax/ndis;
%                 H2(m,n)=H2(m,n) + phixderiv(i+1+(ngrid-abs(r))*ndis)*mu*phixderiv(i+1+ngrid*ndis)*deltax/ndis;
%                 H(m,n)=(k^2)*H1(m,n)+H2(m,n);
%                 
%                 end
%                 
%             elseif n < ngrid+1;                            %unknown boundary elements
%         
%                 for i= MN1(m,n)*ndis:ngrid*ndis-1;         %integration range
%             
%                 T(m,n)=T(m,n) + phix(i+(1+d)+(n-1)*ndis)*rho*phix(i+1+(m-1)*ndis)*deltax/ndis;
%                 H1(m,n)=H1(m,n) + phix(i+(1+d)+(n-1)*ndis)*mu*phix(i+1+(m-1)*ndis)*deltax/ndis;
%                 H2(m,n)=H2(m,n) + phixderiv(i+(1+d)+(n-1)*ndis)*mu*phixderiv(i+1+(m-1)*ndis)*deltax/ndis;    
%                 H(m,n)=(k^2)*H1(m,n)+H2(m,n);
%                 %H2 not completely accurate yet. Although corrected by (1+d)
%                 %contains deviation from correct answer by mag. of a thousandth
%                 %For spline deviation is 0.02-> Because no derivation in 0?
%                 
%                 
%                 
%                 end
%             end
%         end
%     end
% end
% 
% %Last BR nodes:
% for m=IX-ngrid:IX-1;
%     
%     xm=deltax*(m-1);
%     x=xm+(i/ndis)*deltax;
%     rho=a_rho*x+b_rho;
%     mu=a_mu*x+b_mu;
%     
%     for r=-ngrid:IX-m-1;               %range of nodes corrected for boundary
%         
%         for n=m+r;                  %computation of row elements
%         
%             T(m,n)=0.;
%             H1(m,n)=0.;
%             H2(m,n)=0.;
%             H(m,n)=0.;
%             
%             A=zeros(IX-4,IX-4);                 %Arbitrary zero matrix
%             MN1 = [-2 -1 0;-1 -1 0;0 0 0];
%             MN2=blkdiag(A,MN1);                 %matrix describing integration range along upper boundary
%             
%                 if m==n;                        %Kronecker Delta
%                    d=1;
%                 elseif m~=n;
%                        d=0;
%                 end
%                 
%             if n <= IX-ngrid-1;
%             
%                 for i=-(ngrid-abs(r))*ndis:ngrid*ndis-1;      %integration range
%             
%                 T(m,n)=T(m,n) + phix(i+1+(ngrid-abs(r))*ndis)*rho*phix(i+1+ngrid*ndis)*deltax/ndis;
%                 H1(m,n)=H1(m,n) + phix(i+1+(ngrid-abs(r))*ndis)*mu*phix(i+1+ngrid*ndis)*deltax/ndis;
%                 H2(m,n)=H2(m,n) + phixderiv(i+1+(ngrid-abs(r))*ndis)*mu*phixderiv(i+1+ngrid*ndis)*deltax/ndis;
%                 H(m,n)=(k^2)*H1(m,n)+H2(m,n);
%                 
%                 end
%                 
%             elseif n > IX-ngrid-1;
%         
%                 for i= MN2(m,n)*ndis:ngrid*ndis-1;      %integration range
%             
%                 T(m,n)=T(m,n) + phix(i+(1+d)+(IX-n-1)*ndis)*rho*phix(i+1+(IX-m-1)*ndis)*deltax/ndis; 
%                 H1(m,n)=H1(m,n) + phix(i+(1+d)+(IX-n-1)*ndis)*mu*phix(i+1+(IX-m-1)*ndis)*deltax/ndis;
%                 H2(m,n)=H2(m,n) + phixderiv(i+(1+d)+(IX-n-1)*ndis)*mu*phixderiv(i+1+(IX-m-1)*ndis)*deltax/ndis;
%                 H(m,n)=(k^2)*H1(m,n)+H2(m,n);
%                 %H2 not completely accurate yet. Although corrected by (1+d)
%                 %contains deviation from correct answer by mag. of a thousandth
%                 end
%             end
%         end
%     end
% end
        

%_____________________________DEFINITION OF 1D WAVEFORM________________________________

         

% for iomega=iomegamin+1:iomegamax;
%     
%     F(iomega)=((2*iomega^2)/sqrt(pi*iomega^3))*exp((-iomega^2)/iomega^2)/rho;    %rho needs to be position dependant for heterogeneous
%     g(isource,iomega)=F(iomega);
%     
% %     if iomega == 0     % for loop doesnt accept starting at 0
% %        
% %        c(iomega)=0;
% %        
% %     else
%        for m=1:NX-1;
%            for n=1:NX-1;
%             B(m,n,iomega)=((iomega^2)*T(m,n)-H(m,n));
%             %c(iomega)=dot(inv(B(:,:,iomega)),g(:,iomega));      %B is singular...
%           
%             
%            end
%        end
%        
% end
% 
% 
% %LU decomposition
% for iomega=iomegamin+1:iomegamax;
%      d(iomega)=B(m,n,iomega)\g(iomega);
% %     d=B\g;
% end
% 
% for iomega=iomegamin+1:iomegamax;
%     Binv(m,n,iomega) = inv(B(m,n,iomega));
%     d(iomega) = Binv(m,n,iomega)* g(iomega);
% end
% 
% [L,U]=lu(B)
% 
% for iomega=iomegamin+1:iomegamax;
%     d(iomega)=U\(L\g(iomega));
% end

% T=0;
% 





