clc;
clear all;
close all;
warning('off');

%%% Gather Data

TrainFolder_name = 'faces';

Folder_inTS = dir(TrainFolder_name);
no_folders = length(Folder_inTS);

File_Count = 1;
User_Count = 1;

WAITBAR = waitbar(0,'Loading Test Images,Please wait...');

No_Imgs_User = zeros(1,no_folders-2);
for fold_ind = 3:no_folders
    
    waitbar((fold_ind-2)/(no_folders-2))
    
    User_Folder = [TrainFolder_name '\' Folder_inTS(fold_ind).name,'\'];
    ImgList_User = dir(User_Folder);
    No_Imgs_User(User_Count) = length(ImgList_User)-2;
    
    for img_ind = 3:No_Imgs_User(User_Count)+2
        
        Image_Path = User_Folder;
        Image_Name = ImgList_User(img_ind).name;
        Image_Path_Name = [Image_Path,Image_Name];
        
        if strcmp(Image_Name,'Thumbs.db')
            break;
        end
        
        test = imread(Image_Path_Name);
        if length(size(test))==3
            img = imresize(rgb2gray(test),[128,128]);
        else
            img = imresize(test,[128,128]);
        end
        Features = double(imresize(img,[32 32]));
        Data_Matrix(:,File_Count) = Features(:);
        File_Count = File_Count+1;
        
    end
    
    User_Count = User_Count+1;
    
end
close(WAITBAR);

A = Data_Matrix;
Train_Data = A/(diag(sqrt(diag(A'*A))));

%%% Training

Yrec       = sum(A,2)/sum(sqrt(diag(A'*A)));

[Test_File Test_File_Path] = uigetfile('*.jpg;*.pgm;*.png;*.tif','Select a Test Image');
test_image_path = [Test_File_Path Test_File];

test = imread(test_image_path);
if length(size(test))==3
    Test_Image = rgb2gray(test);
else
    Test_Image = test;
end

Features_Test = double(imresize(Test_Image,[32 32]));
Ft_Align = Features_Test(:);
n        = size( Train_Data ,2);

func    = ones(2*n,1);
Aeq     = [Train_Data -Train_Data];
lb      = zeros(2*n,1);
eta     = 0.1;
alpha   = 1;
alphat  = 1;

ln_rat  = 0.5;
max_it  = 10;
Weights = diag(ones( [ 1 , length(Ft_Align) ] ) ); 

WAITBAR = waitbar(0,'Training Test Images,Please wait...');

for iter = 1:max_it
    waitbar(iter/max_it);
    err        = Ft_Align(:).'/sum(sqrt(diag(A'*A))) - Yrec(:).';
    theta      = sum(err)/length(err);
    ln_rat     = ln_rat - 8/theta;
    
    min_sol    = linprog(func,[],[],Aeq,Ft_Align,lb,[],[],[]);
    min_sol    = min_sol(1:n)-min_sol(n+1:2*n);
    
    CS_No_Imgs = cumsum(No_Imgs_User);
    
    tmp        = zeros([1,User_Count-1]);
    for i = 1:User_Count-1
        delta_xi = zeros(length(min_sol),1);
        if i == 1
            delta_xi(1:CS_No_Imgs(i)) = min_sol(1:CS_No_Imgs(i));
        else
            begs = CS_No_Imgs(i-1)+1;
            ends = CS_No_Imgs(i);
            delta_xi(begs:ends) = min_sol(begs:ends);
        end
        tmp(i) = norm(Ft_Align-Train_Data*delta_xi,2);
    end
    
    we_thetha  = exp( ( ln_rat * theta ) - ( ln_rat * err.^2 ) ) ./ ( 1 + exp( ( ln_rat * theta ) - ( ln_rat * err.^2 ) ) );
    Weights    = Weights + diag(we_thetha);    
    alphat     = sum( Weights*( sum(size(repmat(Ft_Align,1,size(Data_Matrix,2))),2) - alpha*sum(Data_Matrix,2) ) ) / ( size(Data_Matrix,1) * size(Data_Matrix,2) );
    
    if(iter==1)
        alpha  = alphat;
    else
        alpha  = alpha + eta*( alphat - alpha );
    end
    
    Yrec       = alpha*Data_Matrix(:,1);
    
end
close(WAITBAR);
clss = find(tmp==min(tmp));
f_name = dir(TrainFolder_name);

msgbox(f_name(clss+2).name);
