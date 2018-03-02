function [L,S] = lps_ist_2D_LLR_plus_S_PCA_V2(param)
%
% L+S reconstruction of undersampled dynamic MRI data using iterative
% soft-thresholding of singular values of L and iterative clipping of
% entries of S (temporal total variation on S)
% 
% Global low rank for first 10 iterations and local low rank afterwards
% Global sparse update throughout
%
% [L,S]=lps_ist_2D_LLR_plus_S(param)
%
% Input variables and reconstruction parameters must be stored in the
% struct param
%
% param.d: undersampled k-t data (nx,ny,nt)
% param.E: data acquisition operator
% param.lambda_L: nuclear-norm weight
% param.lambda_S: temporal weight
%
% Ricardo Otazo (2013) L+S initial code from http://cai2r.net/resources/software/ls-reconstruction-matlab-code
% Tao Zhang (2015) Local Low Rank portion adapted from http://mrsrl.stanford.edu/~tao/software.html
%
% Eugene Milshteyn 01/06/2016
% Update (06/21/2016) global low rank and sparse PCA; was tempfft in first
% 10 iterations

M=param.E'*param.d;
[nx,ny,nt]=size(M);
M=reshape(M,[nx*ny,nt]);
Lpre = M; S=zeros(nx*ny,nt);
ite=0;
bsize = param.bsize;
zpadx = floor(nx/bsize(1));
zpady = floor(ny/bsize(2));
% num_block = [zpadx/bsize(1), zpady/bsize(2)];
z = zeros(bsize(1)*bsize(2),nt-1);
p = zeros(nx*ny,nt-1);


fprintf('\n ********** LLR+S reconstruction **********\n')
% iterations
while(1),
	ite=ite+1;
    if ite<=10
	% low-rank update
        M0=M;
        [Ut,St,Vt]=svd(M-S,0); 
        St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
        L=Ut*St*Vt';
        % sparse update
        [coeff] = pca(M-Lpre);
        S=reshape((SoftThresh((M-Lpre)*coeff,param.lambda_S))*coeff',[nx*ny,nt]);
%         S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[nx,ny,nt]),param.lambda_S)),[nx*ny,nt]);
        % data consistency
        resk=param.E*reshape(L+S,[nx,ny,nt])-param.d;
        M=L+S-reshape(param.E'*resk,[nx*ny,nt]);
        % L_{k-1} for the next iteration
        Lpre=L;
        % print cost function and solution update
        tmp2=param.T*reshape(S,[nx,ny,nt]);
    else
        if ite<=20
            M0=M;
            M = reshape(M,[nx,ny,nt]);
            S = reshape(S,[nx,ny,nt]);
            L = reshape(L,[nx,ny,nt]);
            Y = zeros(zpadx,zpady,nt);
            Z = zeros(zpadx,zpady,nt);
            Y(1:nx, 1:ny, :) = M-S;
            M = reshape(M,[nx*ny,nt]);
            S = reshape(S,[nx*ny,nt]);
    %         Z(1:nx, 1:ny, :) = S;
            randX = randperm(zpadx);
            randY = randperm(zpady);
            Y = circshift(Y, [randX(1)-1, randY(1)-1,  0]);
%             Z = circshift(Z, [randX(1)-1, randY(1)-1,  0]);

            for blockx = 1:zpadx
                    for blocky = 1:zpady

                        tempb = Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:);
                        tempb = reshape(tempb,[bsize(1)*bsize(2),nt]);
                        [Ut,St,Vt] = svd(tempb,'econ');
                        sdiag=diag(SoftThresh(diag(St),St(1)*(param.lambda_L*param.mu)));
                        tempb = Ut*(sdiag)*Vt';
                        tempb = reshape(tempb,[bsize(1),bsize(2),nt]);
                        Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:) = tempb;

                    end
            end

            Y = circshift(Y, [-randX(1)+1, -randY(1)+1,  0]);
    %         Z = circshift(Z, [-randX(1)+1, -randY(1)+1,  0]);
        else
            M0=M;
            M = reshape(M,[nx,ny,nt]);
            S = reshape(S,[nx,ny,nt]);
            L = reshape(L,[nx,ny,nt]);
            Y = zeros(zpadx,zpady,nt);
            Z = zeros(zpadx,zpady,nt);
            Y(1:nx, 1:ny, :) = M-S;
            M = reshape(M,[nx*ny,nt]);
            S = reshape(S,[nx*ny,nt]);
    %         Z(1:nx, 1:ny, :) = S;
            randX = randperm(zpadx);
            randY = randperm(zpady);
            Y = circshift(Y, [randX(1)-1, randY(1)-1,  0]);
    %         Z = circshift(Z, [randX(1)-1, randY(1)-1,  0]);

            for blockx = 1:zpadx
                    for blocky = 1:zpady

                        tempb = Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:);
                        tempb = reshape(tempb,[bsize(1)*bsize(2),nt]);
                        [Ut,St,Vt] = svd(tempb,'econ');
                        sdiag=diag(SoftThresh(diag(St),St(1)*(param.lambda_L*param.mu*0.5)));
                        tempb = Ut*(sdiag)*Vt';
                        tempb = reshape(tempb,[bsize(1),bsize(2),nt]);
                        Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:) = tempb;

                    end
            end

            Y = circshift(Y, [-randX(1)+1, -randY(1)+1,  0]);
    %         Z = circshift(Z, [-randX(1)+1, -randY(1)+1,  0]);
        end
        L = reshape(Y,[nx*ny,nt]);
        % sparse update
        [coeff] = pca(M-Lpre);
        S=reshape((SoftThresh((M-Lpre)*coeff,param.lambda_S))*coeff',[nx*ny,nt]);
%         S_inter=(L*coeff);
%         S = S_inter(:,1:param.lambda_S)*coeff(:,1:param.lambda_S)';
        % data consistency
        resk=param.E*reshape(L+S,[nx,ny,nt])-param.d;
        M=L+S-reshape(param.E'*resk,[nx*ny,nt]);
        % L_{k-1} for the next iteration
        Lpre=L;
        % print cost function and solution update
        tmp2=param.T*reshape(S,[nx,ny,nt]);
    end
    [U,V,W] = svd(L,'econ');

    % print cost function and solution update
%     fprintf('ite: %d, update: %f3\n',ite,norm(M(:)-M0(:))/norm(M0(:)));
    fprintf(' ite: %d , cost: %f3, update: %f3\n', ite,norm(resk(:),2)^2+param.lambda_L*sum(diag(V))+param.lambda_S*norm(diff(S,1,2),1),norm(M(:)-M0(:))/norm(M0(:))); 
    % stopping criteria 
    if (ite > param.nite) || ((norm(M(:)-M0(:))<param.tol*0.1*norm(M0(:))) && ite<11) || (norm(M(:)-M0(:))<param.tol*norm(M0(:)) && ite>=11), break;end
end
L=reshape(L,nx,ny,nt);
S=reshape(S,nx,ny,nt);
end
% soft-thresholding function
function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end    
