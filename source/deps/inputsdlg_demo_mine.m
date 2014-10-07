%INPUTSDLG DEMO (Enhanced input dialog box with multiple data types)

% Written by: Takeshi Ikuma
% Last Updated: May 5 2010
%
% Updated for additional functions F. Hatz 2013

clear; close all;

Title = 'INPUTSDLG Demo Dialog';

%%%% SETTING DIALOG OPTIONS
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'on';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 4; % Horizontal dimension in fields

Prompt = {};
Formats = {};
DefAns = struct([]);

Prompt(1,:) = {['This demo illustrates every type of control that can be placed by INPUTSDLG function ' ... 
   'and demonstrates how Formats input can be used to layout these controls.'],[],[]};
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
Formats(1,1).span = [1 2]; % item is 1 field x 4 fields

Prompt(2,:) = {'Bidder''s Name', 'Name',[]};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'text';
Formats(2,1).size = 200; % automatically assign the height
DefAns(1).Name = 'John Smith';

Prompt(3,:) = {'Bidder''s SSN', 'SSN','(no space or hyphen)'};
Formats(2,2).type = 'edit';
Formats(2,2).format = 'integer';
Formats(2,2).limits = [0 999999999]; % 9-digits (positive #)
Formats(2,2).size = 80;
Formats(2,2).unitsloc = 'bottomleft';
DefAns.SSN = 123456789;

Prompt(4,:) = {'Bidding Price: $', 'Price',[]};
Formats(2,3).type = 'edit';
Formats(2,3).format = 'float';
Formats(2,3).margin = [0 0];
Formats(2,3).limits = [0 inf]; % non-negative decimal number
Formats(2,3).callback = @(hObj,~,~,~)set(hObj,'String',num2str(round(str2double(get(hObj,'String'))*100)/100));
DefAns.Price = 99.99;

Prompt(5,:) = {'Enable enhanced mode' 'EnableEnhancedMode',[]};
Formats(2,4).type = 'check';
DefAns.EnableEnhancedMode = true;

Prompt(end+1,:) = {'Bidder''s Bio File','BioFile',[]};
Formats(3,1).type = 'edit';
Formats(3,1).format = 'file';
Formats(3,1).items = {'*.bio','Biography File (*.bio)';'*.*','All Files'};
Formats(3,1).limits = [0 1]; % single file get
Formats(3,1).size = [-1 0];
Formats(3,1).span = [1 3];  % item is 1 field x 3 fields
d = dir;
files = strcat([pwd filesep],{d(~[d.isdir]).name});
DefAns.BioFile = files{1};


Prompt(end+1,:) = {'Bidder''s Data Folder','DataFolder',[]};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'dir';
Formats(4,1).size = [-1 0];
Formats(4,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.DataFolder = pwd;

Prompt(end+1,:) = {'Save Bidding History To','SaveFile',[]};
Formats(5,1).type = 'edit';
Formats(5,1).format = 'file';
Formats(5,1).limits = [1 0]; % use uiputfile
Formats(5,1).size = [-1 0];
Formats(5,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.SaveFile = files{2};

Prompt(end+1,:) = {'Choose Currency','MoneyUnit',[]};
Formats(7,1).type = 'list';
Formats(7,1).format = 'text';
Formats(7,1).style = 'radiobutton';
Formats(7,1).items = {'U.S. Dollar' 'Euro';'Japanese Yen' ''};
% Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
DefAns.MoneyUnit = 'Japanese Yen';%3; % yen




[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options)
