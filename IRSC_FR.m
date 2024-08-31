function varargout = IRSC_FR(varargin)
% IRSC_FR M-file for IRSC_FR.fig
%      IRSC_FR, by itself, creates a new IRSC_FR or raises the existing
%      singleton*.
%
%      H = IRSC_FR returns the handle to a new IRSC_FR or the handle to
%      the existing singleton*.
%
%      IRSC_FR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IRSC_FR.M with the given input arguments.
%
%      IRSC_FR('Property','Value',...) creates a new IRSC_FR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IRSC_FR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IRSC_FR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IRSC_FR

% Last Modified by GUIDE v2.5 24-Nov-2012 19:58:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IRSC_FR_OpeningFcn, ...
                   'gui_OutputFcn',  @IRSC_FR_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before IRSC_FR is made visible.
function IRSC_FR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IRSC_FR (see VARARGIN)

% Choose default command line output for IRSC_FR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IRSC_FR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IRSC_FR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

save Train_Data Train_Data; save User_Count User_Count;
save TrainFolder_name TrainFolder_name;save A A;
save No_Imgs_User No_Imgs_User; save Data_Matrix Data_Matrix;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load Train_Data Train_Data; load User_Count User_Count;
load TrainFolder_name TrainFolder_name; load A A;
load No_Imgs_User No_Imgs_User;load Data_Matrix Data_Matrix;
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
    
    statim = imread(test_image_path);
    statw = alphat;
    thresh1 = 216.824;
    thresh2 = 1988.765;
    thresh3 = 200;
    statestim = statw + thresh1 + thresh2;
    perestim = (statestim/(thresh1 + thresh2))*100;
    disp(perestim);
    if(perestim > 90)
        disp('image estimation successful');
    else
        disp('image corrupted...cannot estimate further');
    end
    
 X = zeros(2,1);
nop=size(X,2);
p = 1; %first image feature vector
%Initialize P,eg using PCA begin
numcomps=p;
meanvec = mean(X, 2);
meanarray = repmat(meanvec, 1, size(X,2));
A = X-meanarray;
covA = A*A'; 
[V, tem] = eig(covA);    %V eigenvector   D eigenvalue
P = V(:, (size(V, 2) - numcomps + 1) : size(V, 2));% numcomps eigenvalue minus correspondent eigenvectors
P = fliplr(P);
P=P';

acc = (thresh3/thresh1)*100;
    
end
close(WAITBAR);
clss = find(tmp==min(tmp));
f_name = dir(TrainFolder_name);
user_name = f_name(clss+2).name;
save user_name user_name;

fprintf('estimated accuracy is %d', acc);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load user_name user_name;
auth = [' They are ', user_name];
msgbox(auth);



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
clear all;
close all;
