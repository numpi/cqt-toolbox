function cqtguide
%CQTGUIDE Open the online guide for the CQT toolbox.

global cqt_html_file;

if isempty(cqt_html_file)
    fprintf('Compiling the online guide, please wait ... ');
    doc_folder = strcat(fileparts(mfilename('fullpath')),  '/docs');
    addpath(doc_folder);
    cqt_html_file = CompileGuide;
    rmpath(doc_folder);
    fprintf(' done \n');
end

web(cqt_html_file);

