@echo off
set DOCUMENTS_NAMES=main
set TO_DELETE_EXT=-blx.aux -blx.bib .acn .acr .alg .aux .bbl .bcf .blg .cb .cb2 .dvi .fdb_latexmk .fls .fmt .fot .glg .glo .gls .glsdefs .idx .ilg .ind .ist .lof .log .lol .lot .nav .out .pdf .pdfsync .pre .run.xml .snm .sta .synctex .synctex.gz .toc .vrb .xdv
echo.
for %%D in  (%DOCUMENTS_NAMES%) do (
	echo Clean of %%D:
	for %%E in  (%TO_DELETE_EXT%) do (
		if exist %%D%%E (
			echo ^> Removing file %%D%%E
			del %%D%%E >NUL 2>&1
		)
	)
	echo.
)
