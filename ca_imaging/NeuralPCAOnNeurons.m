clear all

addpath(genpath('../ca_source_extraction/utilities'));

nam = 'demoMovie.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels
Yr=reshape(Y,T,d);
%%
close all
rf_x=15;
rf_y=15;
S = 5;
qx=(d1-rf_x)/S+1;
qy=(d2-rf_y)/S+1;
q=qx*qy;
Wmask=zeros(d1,d2,q);
W=zeros(q,d);
counter=1;
for xx=ceil(rf_x/2):S:(d1-floor(rf_x/2))    
    for yy=ceil(rf_y/2):S:(d2-floor(rf_y/2))        
         Wmask(xx-floor(rf_x/2):xx+floor(rf_x/2),yy-floor(rf_y/2):yy+floor(rf_y/2),counter)=rand(rf_x,rf_y);         
         W(counter,:)=reshape(Wmask(:,:,counter),1,[]);
         imagesc(Wmask(:,:,counter))
         disp(counter)
         counter=counter+1;
         drawnow
    end
end
%%
close
M = zeros(q,q);
Y = zeros(q,1);
Ysq=1*10*ones(size(W,1),1);

outer_iter=1;
times_=[];
errors=[];
options.q=q;
for outit=1:outer_iter
    for i = 1:T
        disp(i)
        idx=(outit-1)*T + i;
        tic
        options.gamma=1./Ysq;
        [M,W,Y] = H_AH_NN_PCAFast(M,W,Y, Yr(i,:),options);        
        Ysq = Ysq + Y.^2;
        imagesc(W)                
        drawnow
        %times_(idx,ll)=toc;
%         if mod(idx,nstep_skip_EIGV_errors) ==  0 || idx==n*outer_iter || i == round(n*outer_iter/2)
%             if 0
%                 %disp(['iteration:' num2str(i) ', computing vectors...'])
%                 vectors = orth((pinv(diag(ones(q,1))+M(1:q,1:q))*W(1:q,:))');
%             end
%             %errors(idx,ll)=compute_reconstruction_error(eig_vect,vectors);
%         end
    end
    %             plot(errors(:,ll))
    %             drawnow
end