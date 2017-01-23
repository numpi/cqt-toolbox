function html_file = CompileGuide
%COMPILEGUIDE Compile the toolbox guide from .m files.

doc_files = { ...
	'ArithmeticExample.m' ...
	};

main_file = 'Guide.m';

for i = 1 : length(doc_files)
	publish(doc_files{i});
end

html_file = publish(main_file);


end

