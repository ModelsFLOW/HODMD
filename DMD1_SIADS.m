
 % This file is part of HODMD
 %
 % Copyright (c) 2017 S Le Clainche & J M Vega
 % All rights reserved.
 %
 % Redistribution and use in source and binary forms, with or without
 % modification, are permitted provided that the following conditions
 % are met:
 % 1. Redistributions of source code must retain the above copyright
 %    notice, this list of conditions and the following disclaimer.
 % 2. Redistributions in binary form must reproduce the above copyright
 %    notice, this list of conditions and the following disclaimer in the
 %    documentation and/or other materials provided with the distribution.
 %
 % THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 % ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 % ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 % OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 % HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 % LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 % OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 % SUCH DAMAGE.
 %
 % $FreeBSD$

function  [Vreconst,deltas,omegas,amplitude] =DMD1(V,Time,varepsilon1,varepsilon)

%%%%%%%%%%%%%%%%%%%%%%%%%  DMD-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function solves the DMD-1 algorithm presented in               %%%
%%% Le Clainche & Vega, SIAM J. on Appl. Dyn. Sys. 16(2):882-925, 2017  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %% INPUT: %%
%%% V: snapshot matrix
%%% Time: vector time
%%% varepsilon1: first tolerance (SVD)
%%% varepsilon: second tolerance (DMD-d modes)
%%% %% OUTPUT: %%
%%% Vreconst: reconstruction of the snapshot matrix V
%%% deltas: growht rate of DMD modes
%%% omegas: frequency of DMD modes(angular frequency)
%%% amplitude: amplitude of DMD modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[J,K]=size(V);
[U,Sigma,T]=svd(V,'econ');
sigmas=diag(Sigma);
Deltat=Time(2)-Time(1);
n=length(sigmas);

NormS=norm(sigmas,2);
kk=0;
for k=1:n
    if norm(sigmas(k:n),2)/NormS>varepsilon1
        kk=kk+1;
    end
end
%% Spatial complexity: kk
('Spatial complexity')
kk

U=U(:,1:kk);
%% Create reduced snapshots matrix
hatT=Sigma(1:kk,1:kk)*T(:,1:kk)';
[N,K]=size(hatT);
[hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ');

%% Calculate Koopman operator
hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1';
[Q,MM]=eig(hatR);

eigenvalues=diag(MM);

M=length(eigenvalues);
qq=log(eigenvalues);
deltas=real(qq)/Deltat;
omegas=imag(qq)/Deltat;

%% Calculate amplitudes
Mm=zeros(M*K,M);
Bb=zeros(M*K,1);
for k=1:K
    Mm(1+(k-1)*M:k*M,:)=Q*(MM^(k-1));
    Bb(1+(k-1)*M:k*M,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(Mm,'econ');
a=Vr*(Sigmar\(Ur'*Bb));

u=zeros(M,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end

amplitude=zeros(M,1);
for m=1:M
    aca=U*u(:,m);
    amplitude(m)=norm(aca(:),2)/sqrt(J);
end

UU=[u;deltas';omegas';amplitude']';
UU1=sortrows(UU,-(M+3));

UU=UU1';
u=UU(1:M,:);
deltas=UU(M+1,:);
omegas=UU(M+2,:);
amplitude=UU(M+3,:);
kk2=0;

for m=1:M
    if amplitude(m)/amplitude(1)>varepsilon
        kk2=kk2+1;
    else
    end
end
%% Spectral complexity: number of DMD modes.
('Spectral complexity')
kk2
u=u(:,1:kk2);
deltas=deltas(1:kk2);
omegas=omegas(1:kk2);
amplitude=amplitude(1:kk2);
('Mode number, delta, omega, amplitude')
DeltasOmegAmpl=[deltas',omegas',amplitude']

hatTreconst=zeros(N,K);
for k=1:K
    hatTreconst(:,k)= ContReconst_SIADS(Time(k),Time(1),u,deltas,omegas);
end

Vreconst=U*hatTreconst;






    

