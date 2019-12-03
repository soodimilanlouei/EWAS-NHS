

function summarypanel(X1, Y1, S1, C1)
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c

% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.0625 0.0655487649809055 0.903225814142535 0.910569113444507]);
hold(axes1,'on');

% Create scatter
scatter(X1,Y1,S1,C1,'MarkerFaceAlpha',0.5,'MarkerFaceColor','flat',...
    'MarkerEdgeAlpha',0.5,...
    'MarkerEdgeColor','none');

% Create text
text('Parent',axes1,'FontSize',14,'String','Discretionary liquid fat',...
    'Position',[-0.143682341672302 -0.144092011622887 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Discretionary solid fat',...
    'Position',[0.167533835676719 0.160166607633034 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Animal fat',...
    'Position',[0.155105725638492 0.147169907846187 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Apple juice or cider',...
    'Position',[0.0273181760568339 0.0694908526916057 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Alcohol',...
    'Position',[-0.138336467 -0.115132527 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Synthetic vitamin B6',...
    'Position',[-0.0935357309844844 -0.0140461112727157 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,...
    'String','Added bran from wheat, rice, oat, corn',...
    'Position',[-0.134519336950237 -0.13700437274419 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Natural bran',...
    'Position',[-0.067270249 -0.0769231823386148 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Cereal fiber',...
    'Position',[-0.095601979 -0.0966125788306926 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Cold breakfast cereal',...
    'Position',[-0.0647995666547568 -0.0688767702164328 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Cholesterol',...
    'Position',[0.0931208550819246 0.0936201432457738 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Heme iron',...
    'Position',[0.0585474456338618 0.0771807058800289 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Doughnuts',...
    'Position',[0.078006366 0.0654436901693074 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Myristic acid',...
    'Position',[0.0932677306228869 0.0874706463771131 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Palmitic acid',...
    'Position',[0.140216654828403 0.130195344084214 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Palmitoleic acid',...
    'Position',[0.129233516502135 0.156692816357918 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Stearic acid',...
    'Position',[0.121173106857874 0.135293606930162 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Dietary folate',...
    'Position',[-0.141204944003169 -0.094026435 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,...
    'String',{'Supplemental or fortified','folic acid'},...
    'Position',[-0.137368543313264 -0.0683083688306926 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Natural germ',...
    'Position',[-0.046673089 -0.0897808179904935 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Hamburger',...
    'Position',[0.0807488686445384 0.0333253033957771 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Hotdog',...
    'Position',[0.0711651486082334 0.0620224671194717 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Hydroxyproline',...
    'Position',[0.108181183102587 0.110550050029341 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Liquor',...
    'Position',[-0.0970046959904934 -0.0555888819904934 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,...
    'String','Carbohydrate from milled wholegrain',...
    'Position',[-0.0604324419774896 -0.127297926526347 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Dietary manganese',...
    'Position',[-0.143132245 -0.084703003 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Total manganese',...
    'Position',[-0.137334075479402 -0.077417044 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Animal MUFA',...
    'Position',[0.158882084680159 0.152406514007748 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Plant MUFA',...
    'Position',[-0.118465412 -0.106727055 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Sodium',...
    'Position',[0.135187117622887 0.122558727459038 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Phytate',...
    'Position',[-0.124001832965143 -0.127366420670892 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Peanuts',...
    'Position',[-0.0899600241629697 -0.043062044 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','All processed meats',...
    'Position',[0.0714287408767572 0.0424784872612689 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Raisins or grapes',...
    'Position',[-0.103188203191468 -0.0101108675408528 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Raw carrots',...
    'Position',[-0.0541255808370304 -0.10180093716297 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Red wine',...
    'Position',[-0.108198027696243 -0.081950130156632 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Salad/oil and vinegar dressing',...
    'Position',[-0.167073830146447 -0.0901970997990623 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Total saturated fat',...
    'Position',[0.147263353433112 0.139165979588967 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Supplemental selenium',...
    'Position',[-0.0706595453386148 -0.0808879438370303 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Beverages with sugar',...
    'Position',[0.0805438767264629 0.0379393736978338 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Carbohydrate from wholegrain',...
    'Position',[-0.0542098865636369 -0.105849850237147 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Trans 16:1',...
    'Position',[0.171577107629225 0.168423113736442 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Trans 18:1',...
    'Position',[0.092730772 0.0794729878338615 0],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Trans 18:2',...
    'Position',[0.112137882622887 0.105173048459038 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Alpha tocotrienol',...
    'Position',[-0.0464491578725617 -0.0856755729505852 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Apigenin',...
    'Position',[-0.0966174771701471 -0.04870816961152 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Beta tocotrienol',...
    'Position',[-0.0606605879688111 -0.0948279596717712 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Beta tocopherol',...
    'Position',[-0.0604267665834713 -0.122704936558857 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Isorhamnetin',...
    'Position',[-0.0542218897251967 -0.110245362775819 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Dietary tocopherols',...
    'Position',[-0.129160308606403 -0.131217909880108 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','White bread',...
    'Position',[0.0776286388361508 0.0687565005084525 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','White wine',...
    'Position',[-0.144482646492078 -0.106959153338615 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create text
text('Parent',axes1,'FontSize',14,'String','Yogurt',...
    'Position',[-0.119820518325939 -0.0986587338338615 0],...
    'Color',[0.850980392156863 0.372549019607843 0.00784313725490196]);

% Create ylabel
ylabel('$\beta_{SES}$','FontSize',22,'Interpreter','latex');

% Create xlabel
xlabel('$\beta_{\ no \ SES}$','FontSize',22,'Interpreter','latex');


xlim(axes1,[-0.150000005960464 0.200000002980232]);
ylim(axes1,[-0.150000005960464 0.200000002980232]);

axis(axes1,'square');
% Set the remaining axes properties
set(axes1,'FontSize',20);
% Create line
annotation(figure1,'line',[0.631359181535929 0.60222384319006],...
    [0.574199061468049 0.610784427321717],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.326801144652677 0.320222197284253],...
    [0.286691578317208 0.409659057992008],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.302238805970149 0.297015860546889],...
    [0.269992082343626 0.317995288028889],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.295398009950249 0.285115774510904],...
    [0.266825019794141 0.304441129328158],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.311567164179104 0.306110425317023],...
    [0.269992082343626 0.337067666994961],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.307274716642351 0.322761194029851],...
    [0.252894289631735 0.254948535233571],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.304852616616166 0.317786069651741],...
    [0.245698395226236 0.245447347585115],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.292952530580182 0.259950248756219],...
    [0.254856806287778 0.276326207442597],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.297015860546888 0.364427860696517],...
    [0.237013112411413 0.229612034837688],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create textbox
annotation(figure1,'textbox',...
    [0.181451127819549 0.866869918699187 0.277195488721805 0.0345528455284689],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'String','$\bullet$ Both analyses',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.181451127819548 0.83130081300814 0.277195488721805 0.0345528455284689],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725],...
    'String','$\bullet$ ``No SES'''' only',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create line
annotation(figure1,'line',[0.289786967418546 0.253109452736318],...
    [0.249657223413089 0.2541567695962],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.265547263681592 0.246890547263681],...
    [0.226444972288203 0.231195566112431],...
    'Color',[0.458823529411765 0.43921568627451 0.701960784313725],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.274487524782104 0.350124378109453],...
    [0.214773670895856 0.191607284243864],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.28893128343246 0.347636815920398],...
    [0.234637815499299 0.216152019002375],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.231095462536938 0.245646766169154],...
    [0.214051908927644 0.209817893903405],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

% Create line
annotation(figure1,'line',[0.262835259791269 0.337083380092021],...
    [0.194672962169051 0.141827433713767],...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274],...
    'LineStyle','--');

set(gcf, 'Position', [62, 68, 1597,1263])




function HRSES(data1)
%  DATA1:  histogram data

% Create figure
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create histogram
histogram(data1,'LineWidth',1,...
    'EdgeColor',[0 0.447058823529412 0.741176470588235],...
    'FaceColor','none',...
    'Normalization','pdf',...
    'BinMethod','auto');

% Create ylabel
ylabel('$PDF$','Interpreter','latex');

% Create xlabel
xlabel('$HR_{SES}$','Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure SI4
opts = delimitedTextImportOptions("NumVariables", 13);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["varname", "Exposure", "coef", "exp_coef", "se_coef", "z", "p_val", "LB", "UB",  "PH_pval", "VIF", "FDR"];
opts.VariableTypes = ["char", "char", "double", "double", "double", "double", "double", "double", "double",  "double", "double", "double"];
opts = setvaropts(opts, [1, 2], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
NoSES = readtable("EWAS_results.csv", opts);
% Clear temporary variables
clear opts

opts = delimitedTextImportOptions("NumVariables", 14);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["varname", "Exposure", "coef", "exp_coef", "se_coef", "z", "p_val", "LB", "UB", "PH_pval", "VIF", "type", "FDR"];
opts.VariableTypes = ["char", "char", "double", "double", "double", "double", "double", "double", "double",  "double", "double", "categorical", "double"];
opts = setvaropts(opts, [1, 2], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 13], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
SES = readtable("EWAS_results_withSES.csv", opts);
% Clear temporary variables
clear opts

SES = sortrows(SES,'varname','ascend');
NoSES = sortrows(NoSES,'varname','ascend');

% summary table
summarytable=table;
summarytable.varname=NoSES.varname;
summarytable.Exposure=NoSES.Exposure;
summarytable.type=SES.type;
summarytable.coef_old=NoSES.coef;
summarytable.coef_new=SES.coef;
summarytable.HR_old=NoSES.exp_coef;
summarytable.HR_new=SES.exp_coef;
summarytable.p_val_old=NoSES.p_val;
summarytable.p_val_new=SES.p_val;
summarytable.FDR_old=NoSES.FDR;
summarytable.FDR_new=SES.FDR;
summarytable.filter_old=(summarytable.FDR_old<0.05);
summarytable.filter_new=(summarytable.FDR_new<0.05);

label=cell(height(summarytable),1);
label(:)={'none'};
label((summarytable.filter_old+summarytable.filter_new)==2)={'both'};
label((summarytable.filter_old-summarytable.filter_new)==1)={'no SES'};
label((summarytable.filter_new-summarytable.filter_old)==1)={'SES'};
summarytable.label=label;

cmap=[220,220,220;217,95,2;117,112,179]/255;
label=categorical(label);
cmv=zeros(length(cmv), 3);
cmv(label=='none', :)=repmat(cmap(1,:),sum(label=='none'),1);
cmv(label=='both', :)=repmat(cmap(2,:),sum(label=='both'),1);
cmv(label=='no SES', :)=repmat(cmap(3,:),sum(label=='no SES'),1);

summarypanel(summarytable.coef_old,summarytable.coef_new,[] , cmv)
%% Figure SI5
opts = delimitedTextImportOptions("NumVariables", 9);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["VarName1", "coef", "exp_coef", "se_coef", "z", "p_val", "LB", "UB", "PH_pval"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
EWAS_results_withSES = readtable("EWAS_results_withSES.csv", opts);
% Clear temporary variables
clear opts



HRSES([EWAS_results_withSES.exp_coef])
disp("average HRSES")
disp(mean([EWAS_results_withSES.exp_coef]))
disp("std HRSES")
disp(std([EWAS_results_withSES.exp_coef]))


