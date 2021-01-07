% GUI_STD_ENVTOPO MATLAB code for gui_std_envtopo.fig
%      GUI_STD_ENVTOPO, by itself, creates a new GUI_STD_ENVTOPO or raises the existing
%      singleton*.
%
%      H = GUI_STD_ENVTOPO returns the handle to a new GUI_STD_ENVTOPO or the handle to
%      the existing singleton*.
%
%      GUI_STD_ENVTOPO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_STD_ENVTOPO.M with the given input arguments.
%
%      GUI_STD_ENVTOPO('Property','Value',...) creates a new GUI_STD_ENVTOPO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_std_envtopo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_std_envtopo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_std_envtopo

% Last Modified by GUIDE v2.5 14-Aug-2018 14:38:18

% History
% 06/27/2019 Makoto. Improved efficiency.
% 05/29/2019 Makoto. Returns error when 'Outlier Cls' is detected. Sorry Johanna.
% 09/18/2019 Makoto. Conversion from numeric condition labels to strings supported. 
% 08/14/2019 Makoto. PVAF/PPAF supported.
% 08/10/2019 Makoto. Updated to be compatible with EEGLAB15.
% 03/27/2019 Clement and Makoto. Created.



function varargout = gui_std_envtopo(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_std_envtopo_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_std_envtopo_OutputFcn, ...
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



% --- Executes just before gui_std_envtopo is made visible.
function gui_std_envtopo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_std_envtopo (see VARARGIN)

% Load STUDY and ALLEEG
STUDY  = evalin('base', 'STUDY');
ALLEEG = evalin('base', 'ALLEEG');

% Check if there is 'Outliers Cls...' (05/29/2019 Makoto.)
finalClusterName = STUDY.cluster(end).name;
if any(strfind(finalClusterName, 'Outliers Cls'))
    error('''Outlier Cls XX'' detected, which is not supported: Use ''outlier 2'' instead, by enabling ''Separate outliers (enter std.)''')
end

% Check if EEGLAB15 >=3 variables. 08/10/2018 Makoto
if length(STUDY.design(STUDY.currentdesign).variable)>=3
    error('Three-way interaction not supported. Please use one- or two-way.')
end
    
% Set the default window length to be that of the data. % 08/10/2018 Makoto. Updated to be compatible with EEGLAB15.
if ~isfield(STUDY.cluster, 'erpdata')
    try   % EEGLAB 14
        STUDY = std_readerp(STUDY, ALLEEG, 'design', STUDY.currentdesign, 'clusters', 2:length(STUDY.cluster), 'singletrials', 'off', 'timerange', [ALLEEG(1,1).xmin*1000 ALLEEG(1,1).xmax*1000]);
        erpTimes = STUDY.cluster(1,2).erptimes;
        assignin('base', 'STUDY', STUDY);
    catch % EEGLAB 2019
        disp(sprintf('\n!!! EEGLAB version 15 or later detected !!!\n'))
        [~, ~, erpTimes] = std_readdata(STUDY, ALLEEG, 'design', STUDY.currentdesign, 'clusters', 2, 'datatype', 'erp', 'singletrials', 'off', 'timerange', [ALLEEG(1,1).xmin*1000 ALLEEG(1,1).xmax*1000]);
    end
else     % EEGLAB 14
    erpTimes = STUDY.cluster(1,2).erptimes;
end

% Set strings in GUI.
set(handles.timeRangeToPlot, 'string', sprintf('%.0f %.0f', [min(erpTimes) max(erpTimes)]))
if floor(min(abs(erpTimes))) == 0
    set(handles.timeRangeToRankClusters, 'string', sprintf('%.0f %.0f', [0 max(erpTimes)]))
else
    set(handles.timeRangeToRankClusters, 'string', sprintf('%.0f %.0f', [min(erpTimes) max(erpTimes)]))
end

% Choose default command line output for gui_std_envtopo
handles.output = hObject;

% Populate popup menus with var1 and var2.
var1 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
if iscell([var1{:}]) % Decode combined names by STUDY design (10/20/2016 Makoto)
    for n = 1:length(var1)
        tmpName = var1{n};
        if iscell(tmpName) % This means tmpName contains multiple cells containing names.
            newTmpName = '';
            for m = 1:length(tmpName)
                newTmpName = [newTmpName tmpName{m}];
                if m ~= length(tmpName)
                    newTmpName = [newTmpName '+'];
                end
            end
            newVar1{n} = newTmpName;
        else
            newVar1{n} = tmpName;
        end
    end
    var1 = newVar1;
end


if length(STUDY.design(STUDY.currentdesign).variable) == 1 % 08/10/2018 Makoto. Updated to be compatible with EEGLAB15. This is the case in EEGLAB15 with only one condition.
    var2 = {''};
else
    var2 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
    if iscell([var2{:}]) % Decode combined names by STUDY design (10/20/2016 Makoto)
        for n = 1:length(var2)
            tmpName = var2{n};
            if iscell(tmpName) % This means tmpName contains multiple cells containing names.
                newTmpName = '';
                for m = 1:length(tmpName)
                    newTmpName = [newTmpName tmpName{m}];
                    if m ~= length(tmpName)
                        newTmpName = [newTmpName '+'];
                    end
                end
                newVar2{n} = newTmpName;
            else
                newVar2{n} = tmpName;
            end
        end
        var2 = newVar2;
    end
end

% Convert var1 and var2 to strings. 09/17/2018 Makoto.
if ~isempty(var1)
    for loopIdx = 1:length(var1)
        if isnumeric(var1{loopIdx})
            var1{loopIdx} = num2str(var1{loopIdx});
        end
    end
end
if ~isempty(var2)
    for loopIdx = 1:length(var2)
        if isnumeric(var2{loopIdx})
            var2{loopIdx} = num2str(var2{loopIdx});
        end
    end
end


kk = 0; designs = {};
for ii = 1:length(var1)
    for jj = 1:length(var2)
        kk = kk+1;
        
        % 08/10/2018 Makoto. Small improvement not to show '_' when var2 is empty.
        if isempty(var2{1,1})
            designs(kk) = cellstr(strcat(var1(ii)));
        else
            designs(kk) = cellstr(strcat(var1(ii),'_', var2(jj)));
        end
    end
end
varCombinations = ['Off'; designs'];

set(handles.singleConditionDropdown,'String', varCombinations);
set(handles.subtraction1Dropdown,'String', varCombinations);
set(handles.subtraction2Dropdown,'String', varCombinations);
set(handles.interaction1Dropdown,'String', varCombinations);
set(handles.interaction2Dropdown,'String', varCombinations);
set(handles.interaction3Dropdown,'String', varCombinations);
set(handles.interaction4Dropdown,'String', varCombinations);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_std_envtopo wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = gui_std_envtopo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function timeRangeToPlot_Callback(hObject, eventdata, handles)
% hObject    handle to timeRangeToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeRangeToPlot as text
%        str2double(get(hObject,'String')) returns contents of timeRangeToPlot as a double



% --- Executes during object creation, after setting all properties.
function timeRangeToPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeRangeToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function timeRangeToRankClusters_Callback(hObject, eventdata, handles)
% hObject    handle to timeRangeToRankClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeRangeToRankClusters as text
%        str2double(get(hObject,'String')) returns contents of timeRangeToRankClusters as a double



% --- Executes during object creation, after setting all properties.
function timeRangeToRankClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeRangeToRankClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numberOfLargestContributingClusters_Callback(hObject, eventdata, handles)
% hObject    handle to numberOfLargestContributingClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numberOfLargestContributingClusters as text
%        str2double(get(hObject,'String')) returns contents of numberOfLargestContributingClusters as a double



% --- Executes during object creation, after setting all properties.
function numberOfLargestContributingClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numberOfLargestContributingClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elsePlotTheseClusterNumbers_Callback(hObject, eventdata, handles)
% hObject    handle to elsePlotTheseClusterNumbers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elsePlotTheseClusterNumbers as text
%        str2double(get(hObject,'String')) returns contents of elsePlotTheseClusterNumbers as a double



% --- Executes during object creation, after setting all properties.
function elsePlotTheseClusterNumbers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elsePlotTheseClusterNumbers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clusterNumbersToExclude_Callback(hObject, eventdata, handles)
% hObject    handle to clusterNumbersToExclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusterNumbersToExclude as text
%        str2double(get(hObject,'String')) returns contents of clusterNumbersToExclude as a double



% --- Executes during object creation, after setting all properties.
function clusterNumbersToExclude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterNumbersToExclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function optionalInputs_Callback(hObject, eventdata, handles)
% hObject    handle to optionalInputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of optionalInputs as text
%        str2double(get(hObject,'String')) returns contents of optionalInputs as a double



% --- Executes during object creation, after setting all properties.
function optionalInputs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optionalInputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numIterationEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numIterationEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numIterationEdit as text
%        str2double(get(hObject,'String')) returns contents of numIterationEdit as a double



% --- Executes during object creation, after setting all properties.
function numIterationEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numIterationEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in singleConditionDropdown.
function singleConditionDropdown_Callback(hObject, eventdata, handles)
% hObject    handle to singleConditionDropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns singleConditionDropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from singleConditionDropdown



% --- Executes during object creation, after setting all properties.
function singleConditionDropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleConditionDropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in subtraction1Dropdown.
function subtraction1Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to subtraction1Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subtraction1Dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subtraction1Dropdown
% hObject.String = handles.singleConditionDropdown.String;



% --- Executes during object creation, after setting all properties.
function subtraction1Dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subtraction1Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in subtraction2Dropdown.
function subtraction2Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to subtraction2Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subtraction2Dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subtraction2Dropdown
% hObject.String = handles.singleConditionDropdown.String;



% --- Executes during object creation, after setting all properties.
function subtraction2Dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subtraction2Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in interaction1Dropdown.
function interaction1Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to interaction1Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns interaction1Dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from interaction1Dropdown



% --- Executes during object creation, after setting all properties.
function interaction1Dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interaction1Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in interaction2Dropdown.
function interaction2Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to interaction2Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns interaction2Dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from interaction2Dropdown



% --- Executes during object creation, after setting all properties.
function interaction2Dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interaction2Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in interaction3Dropdown.
function interaction3Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to interaction3Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns interaction3Dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from interaction3Dropdown



% --- Executes during object creation, after setting all properties.
function interaction3Dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interaction3Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in interaction4Dropdown.
function interaction4Dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to interaction4Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns interaction4Dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from interaction4Dropdown



% --- Executes during object creation, after setting all properties.
function interaction4Dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interaction4Dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pvafStatVar.
function pvafStatVar_Callback(hObject, eventdata, handles)
% hObject    handle to pvafStatVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pvafStatVar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pvafStatVar



% --- Executes during object creation, after setting all properties.
function pvafStatVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pvafStatVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pophelp('std_envtopo.m');



% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)

% Generate 'arguments' from Inputs.
STUDY  = evalin('base', 'STUDY');
ALLEEG = evalin('base', 'ALLEEG');

% 08/10/2018 Makoto. Updated to be compatible with EEGLAB15. This is the case in EEGLAB15 with only one condition.
if length(STUDY.design(STUDY.currentdesign).variable) == 1;
    var2 = {''};
else
    var2 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
end

% Prepare options.
options = '';
if ~isempty(get(handles.timeRangeToPlot,'String')),                     options = [ options '''timerange'',[' get(handles.timeRangeToPlot,'string') ']' ]; end;
if ~isempty(get(handles.timeRangeToRankClusters,'String')),             options = [ options ',''limitcontribtime'',[' get(handles.timeRangeToRankClusters,'string') ']' ]; end;
if ~isempty(get(handles.numberOfLargestContributingClusters,'String')), options = [ options ',''clusterIdxToUse'',[-' get(handles.numberOfLargestContributingClusters,'string') ']' ]; end;
if ~isempty(get(handles.elsePlotTheseClusterNumbers,'String')),         options = [ options ',''clusterIdxToUse'',[' get(handles.elsePlotTheseClusterNumbers,'string') ']' ]; end;
if ~isempty(get(handles.clusterNumbersToExclude,'String')),             options = [ options ',''clust_exclude'',[' get(handles.clusterNumbersToExclude,'string') ']' ]; end;
if ~isempty(get(handles.numIterationEdit,'String'));                    options = [ options ',''statPvaf'',[' get(handles.numIterationEdit,'string') ']']; end; 
if ~isempty(get(handles.optionalInputs,'String')),                      options = [ options ',' get(handles.optionalInputs, 'string') ]; end;
switch get(handles.selectMeasurePopupmenu, 'value')
    case 1
        options = [ options ', ''measure'', ''pvaf'''];
    case 2
        options = [ options ', ''measure'', ''ppaf'''];
end


% Define arglist.designVector from dropdown menus
singleCondition  = get(handles.singleConditionDropdown ,'Value');
subtraction1 = get(handles.subtraction1Dropdown,'Value');
subtraction2 = get(handles.subtraction2Dropdown,'Value');
interaction1 = get(handles.interaction1Dropdown,'Value');
interaction2 = get(handles.interaction2Dropdown,'Value');
interaction3 = get(handles.interaction3Dropdown,'Value');
interaction4 = get(handles.interaction4Dropdown,'Value');
argdesignVectorList  = [singleCondition, subtraction1, subtraction2, interaction1, interaction2, interaction3, interaction4];
argdesignVectorListCheck = length(find((argdesignVectorList ~= 1)));
if singleCondition ~= 1 && argdesignVectorListCheck == 1
    designVector(1,1) = ceil((singleCondition - 1)/length(var2));
    designVector(1,2) = mod(singleCondition-1-1, length(var2))+1;
elseif subtraction1 ~= 1 && subtraction2 ~= 1 && argdesignVectorListCheck == 2
    designVector(1,1) = ceil((subtraction1 - 1)/length(var2));
    designVector(1,2) = mod(subtraction1-1-1, length(var2))+1;
    designVector(1,3) = ceil((subtraction2 - 1)/length(var2));
    designVector(1,4) = mod(subtraction2-1-1, length(var2))+1;
elseif interaction1 ~= 1 && interaction2 ~= 1 && ...
        interaction3 ~= 1 && interaction4 ~= 1 && argdesignVectorListCheck == 4
    designVector(1,1) = ceil((interaction1 - 1)/length(var2));
    designVector(1,2) = mod(interaction1-1-1, length(var2))+1;
    designVector(1,3) = ceil((interaction2 - 1)/length(var2));
    designVector(1,4) = mod(interaction2-1-1, length(var2))+1;
    designVector(1,5) = ceil((interaction3 - 1)/length(var2));
    designVector(1,6) = mod(interaction3-1-1, length(var2))+1;
    designVector(1,7) = ceil((interaction4 - 1)/length(var2));
    designVector(1,8) = mod(interaction4-1-1, length(var2))+1;
else
    sprintf('Please select Combination, Subtraction, or Interaction and turn other fields "Off"')
end
options = [ options ',''designVector'',[' num2str(designVector) ']' ];


if length(STUDY.cluster) == 1
    errordlg2('Cannot plot envtopo with the parent cluster only');
end

% Update STUDY in the base workspace
assignin('base', 'STUDY', STUDY);

% Perform std_envtopo
arguments = eval([ '{' options '}' ]);
tic;
STUDY = std_envtopo(STUDY, ALLEEG, arguments{:});
measuredTime = toc;

% Report results
resultReport = sprintf('std_envtopo computation time: %3.0f seconds.',round(measuredTime));
disp(resultReport)


% --- Executes on selection change in selectMeasurePopupmenu.
function selectMeasurePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to selectMeasurePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectMeasurePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectMeasurePopupmenu


% --- Executes during object creation, after setting all properties.
function selectMeasurePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectMeasurePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
