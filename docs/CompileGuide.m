function html_file = CompileGuide
%COMPILEGUIDE Compile the toolbox guide from .m files.

if exist('stylesheet.xsl', 'file')
    publish_cmd = @(file) publish(file, 'stylesheet', 'stylesheet.xsl');
else
    publish_cmd = @(file) publish(file);
end

doc_files = { ...
    'QTDefinition.m', ...
    'Arithmetic.m', ...
    'LinearSystems.m', ...
    'JacksonExample.m', ...
    'SquareRoot.m', ...
    'QuadraticExample.m', ...
    'doc_eig.m' ...
    };

main_file = 'index.m';

for i = 1 : length(doc_files)
    publish_cmd(doc_files{i});
end

html_file = publish_cmd(main_file);

if exist('numpi.css', 'file')
    copyfile('numpi.css', 'html/numpi.css');
end


end

