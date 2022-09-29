%######################################################################
%######################################################################
%IFM Template fuer Uebungsblaetter
%######################################################################
%######################################################################
%
%by Robin Pfefferkorn 27.10.2020
%
%Changelog:
%-tt.mm.jjjj Name of Editor: Changes

%######################################################################
%######################################################################

%-----------------------------------------------------------------------
%preliminaries of new document class
%-----------------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ifmdocument}[2020/10/27 basic ifm class for simple documents]

%very basic packages
\RequirePackage{fix-cm}
\RequirePackage{kvoptions}
\RequirePackage{anyfontsize}

%-----------------------------------------------------------------------
%options of class
%-----------------------------------------------------------------------
\SetupKeyvalOptions{family=ifmdocument,prefix=ifmdocument@}

%language
\DeclareStringOption[german]{language}
\DeclareVoidOption{german}{\ClassError{ifmdocument}
	{set language in document option using 'language=german', not 'german'}}
\DeclareVoidOption{english}{\ClassError{ifmdocument}
	{set language in document option using 'language=english', not 'english'}}
%option for header and footer
\DeclareBoolOption[false]{professorsinheader}
\DeclareBoolOption[false]{pagenumberinfooter}
\DeclareBoolOption[false]{authorinfooter}
\DeclareBoolOption[false]{dateinfooter}
\DeclareBoolOption[false]{sectioninheader}
\DeclareBoolOption[false]{pagenumberinheader}
%biblography
\DeclareStringOption[bibtex]{bibbackend}
\DeclareStringOption[none]{bibstyle}

%all other options to koma-script
\DeclareDefaultOption{\PassOptionsToClass{\CurrentOption}{scrartcl}}

%process options
\ProcessKeyvalOptions{ifmdocument}
\relax

%set language option
%----------------------------
\def\ifmdocument@languagegerman{german}
\def\ifmdocument@languageenglish{english}

\ifx\ifmdocument@language\ifmdocument@languagegerman
	\usepackage[english,ngerman]{babel}
	\AtBeginDocument{\selectlanguage{ngerman}}
\else\ifx\ifmdocument@language\ifmdocument@languageenglish
	\usepackage[ngerman,english]{babel}
	\AtBeginDocument{\selectlanguage{english}}
\else
	\ClassError{ifmdocument}{language unavailable. language must be 'english' or 'german'}{}
\fi\fi

%set bibbackend option
%----------------------------
\def\ifmdocument@bibbackendbiber{biber}
\def\ifmdocument@bibbackendbibtex{bibtex}

\ifx\ifmdocument@bibbackend\ifmdocument@bibbackendbibtex
	\PassOptionsToPackage{backend=bibtex}{biblatex}
\else\ifx\ifmdocument@bibbackend\ifmdocument@bibbackendbiber
	\PassOptionsToPackage{backend=biber}{biblatex}
\else
	\PackageWarning{ifmdocument}{ifmdocument: bibbackend unknown - default ('biber') selected. known input 'bibtex' or 'biber'}
	\PassOptionsToPackage{backend=biber}{biblatex}
\fi\fi

%set bibstyle
%----------------------------
\def\ifmdocument@bibstylenone{none}
\def\ifmdocument@bibstylenumeric{numeric}
\def\ifmdocument@bibstylealpha{alpha}
\def\ifm@bibstyleauthoryear{authoryear}

\ifx\ifmdocument@bibstyle\ifmdocument@bibstylenone
	%no bibliography(standard)
\else\ifx\ifmdocument@bibstyle\ifmdocument@bibstylenumeric
	\PassOptionsToPackage{style=numeric-comp,bibstyle=numeric}{biblatex}
\else\ifx\ifmdocument@bibstyle\ifmdocument@bibstylealpha
	\PassOptionsToPackage{style=alphabetic,bibstyle=alphabetic}{biblatex}
\else\ifx\ifm@bibstyle\ifm@bibstyleauthoryear
\PassOptionsToPackage{style=authoryear,bibstyle=authoryear,uniquename=init,dashed=false}{biblatex}
\else
	\ClassWarning{ifmdocument}{ifmdocument: bibstyle unsupported. may cause errors in biblatex}
	\PassOptionsToPackage{style=\ifmdocument@bibstyle,bibstyle=\ifmdocument@bibstyle}{biblatex}
\fi\fi\fi\fi

%-----------------------------------------------------------------------
%based on KOMAscript article
%-----------------------------------------------------------------------
\LoadClass{scrartcl}

%----------------------------------------------------------------------
%various packages
%----------------------------------------------------------------------
%encoding
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsthm}
\usepackage[utf8]{inputenc}
\usepackage{xcolor}
\usepackage{psfrag}
\usepackage{graphicx}
\usepackage{lmodern}
\usepackage{csquotes}
\usepackage{blindtext}

%output of units
\usepackage{nicefrac}
\usepackage[binary-units=true]{siunitx}
\sisetup{per-mode=fraction,fraction-function=\nicefrac,output-decimal-marker={,}}

%-----------------------------------------------------------------------
%general content (like author etc.)
%-----------------------------------------------------------------------
%fallbacks if not all fields are defined
\gdef\@headerblock{}
\ifx\ifmdocument@language\ifmdocument@languagegerman
	\gdef\@kitlogopath{kitlogo}
\else\ifx\ifmdocument@language\ifmdocument@languageenglish
	\gdef\@kitlogopath{kitlogoEN}
\fi\fi

%macros for definitions
\newcommand{\headerblock}[1]{%
	\ififmdocument@professorsinheader%
		\PackageError{ifmdocument}{if ifmdocument option 'professorinheader' is set, command 'headerblock' may not be used}{}
	\else
		\gdef\@headerblock{#1}
	\fi
}
\newcommand{\kitlogopath}[1]{\gdef\@kitlogopath{#1}}
\AtBeginDocument{
	\newcommand*{\thetitle}{\@title}
	\newcommand*{\theauthor}{\@author}
	\newcommand*{\thedate}{\@date}
}


%-----------------------------------------------------------------------
%Fußzeile und kopfzeile, formatting
%-----------------------------------------------------------------------
\usepackage{scrlayer-scrpage}
\pagestyle{scrheadings}
\clearpairofpagestyles
\raggedbottom
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5em}
%formatieren
\setkomafont{pagefoot}{\normalfont\small}
\setkomafont{pagehead}{\normalfont\scshape\small}
\setkomafont{pagenumber}{\normalfont\small}
\addtokomafont{subsection}{\normalsize}
\addtokomafont{section}{\large}
\RedeclareSectionCommands[afterskip=0.5pt,beforeskip=5pt]{subsection}

%optional settings for footer 
\ififmdocument@pagenumberinfooter
	\ofoot*{\pagemark}
\fi
\ififmdocument@authorinfooter
	\def\@footerauthor{\theauthor}
\else
	\def\@footerauthor{}
\fi
\ififmdocument@dateinfooter
	\def\@footerdatetmp{\the\day.\the\month.\the\year}
	\ifx\@footerauthor\empty
		\def\@footerdate{\@footerdatetmp}
	\else
		\def\@footerdate{, \@footerdatetmp}
	\fi
\else
	\def\@footerdate{}
\fi
\ifoot*{\@footerauthor\@footerdate}

%optional settings for header
\newif\ififm@usehead\ifm@useheadfalse
\ififmdocument@sectioninheader
	\automark[subsection]{section}
	\rehead[]{\leftmark}
	\lohead[]{\rightmark}
	\ifm@useheadtrue
\fi
\ififmdocument@pagenumberinheader
	\lehead[]{\pagemark}
	\rohead[]{\pagemark}
	\ifm@useheadtrue
\fi
\ififm@usehead
	\KOMAoption{headsepline}{0.5pt}
\fi
	
%----------------------------------------------------------------------
%page style
%----------------------------------------------------------------------
%heights for header and footer
\newlength{\ifm@headsep}
\newlength{\ifm@head}
\newlength{\ifm@foot}
\newlength{\ifm@top}

%standard heigths
\setlength{\ifm@headsep}{5mm}
\setlength{\ifm@head}{8mm}
\setlength{\ifm@foot}{11mm}

%if no header is active set lengths to zero
\newcommand{\ifm@inclheadfoot}{%
	\ififm@usehead%
		includeheadfoot%
	\else%
		includefoot%
	\fi%
}

\def\@geomoptions{a4paper,left=20mm,right=20mm,top=20mm,bottom=11mm,foot=\ifm@foot,head=\ifm@head,headsep=\ifm@headsep,\ifm@inclheadfoot}
\usepackage[\@geomoptions]{geometry}
\KOMAoptions{usegeometry}	

%----------------------------------------------------------------------
%Hyperref
%----------------------------------------------------------------------
%set colorlinks to true for versio w/o coloredlinks and border
\usepackage[pdfstartpage=1,
colorlinks=false,citecolor=black,linkcolor=black,urlcolor=black,
citebordercolor=green,linkbordercolor=red,breaklinks=true]{hyperref}
%remove underlines from pdfstrings (pdfauthor,...)
\pdfstringdefDisableCommands{%
	\def\underline#1{<#1>}%
}
\AtBeginDocument{
	\hypersetup{pdftitle={{\@title}},pdfauthor={\@author}}
}

%-----------------------------------------------------------------------
%tikz
%-----------------------------------------------------------------------
\usepackage{tikzscale,pgfplots}
\pgfplotsset{compat=newest}
\usepgfplotslibrary{patchplots}
\newcommand\figurescale{1}
\newlength\figH
\newlength\figW

%----------------------------------------------------------------------
%Literaturverzeichnis
%----------------------------------------------------------------------
\ifx\ifmdocument@bibstyle\ifmdocument@bibstylenone\else
%select biblatex with correct options
\usepackage[natbib=true]{biblatex}
\ExecuteBibliographyOptions{sorting=nyt,maxbibnames=10,maxcitenames=2,uniquename=false}
\setcounter{biburlnumpenalty}{9000}
\setcounter{biburlucpenalty}{9000}

\DeclareNameFormat{author}{%
	\renewcommand*{\multinamedelim}{\addcomma\addspace}%
	\nameparts{#1}%
	\usebibmacro{name:family-given}
	{\textsc{\namepartfamily}}
	{\namepartgiveni}
	{\namepartprefix}
	{\namepartsuffix}%
	\usebibmacro{name:andothers}%
}
\DeclareNameFormat{editor}{%
	\renewcommand*{\multinamedelim}{\addcomma\addspace}%
	\nameparts{#1}%
	\usebibmacro{name:family-given}
	{\textsc{\namepartfamily}}
	{\namepartgiveni}
	{\namepartprefix}
	{\namepartsuffix}%
	\usebibmacro{name:andothers}%
}

\renewbibmacro*{journal+issuetitle}{%
	\usebibmacro{journal}%
	\setunit*{\addcomma\addspace}%
	\iffieldundef{series}
	{}
	{\newunit
		\printfield{series}%
		\setunit{\addspace}}%
	\usebibmacro{volume+number+eid}%
	\setunit{\addspace}%
	\usebibmacro{issue+date}%
	\setunit{\addcolon\space}%
	\usebibmacro{issue}%
	\newunit
}
\renewbibmacro*{volume+number+eid}{%
	\printfield{volume}%
	\iffieldundef{number}{}{
		\unspace\printtext[parens]{\printfield{number}}}%
	\setunit*{\addcolon\addspace}
	\printfield{pages}%
	\addcomma}
\renewbibmacro*{issue+date}{%
	\printtext{%
		\printfield{issue}%
		\setunit*{\addspace}%
		\usebibmacro{date}}%
	\newunit}
\renewbibmacro*{note+pages}{%
	\printfield{note}%
	\newunit}

\DefineBibliographyStrings{ngerman}{
	andothers = {{et\,al\adddot}},
}
\DeclareFieldFormat[article]{pages}{#1\isdot}
\DeclareFieldFormat[article]{title}{\addspace#1\isdot}
\DeclareFieldFormat[incollection]{title}{\addspace#1\isdot}
\DeclareFieldFormat[inbook]{title}{\addspace#1\isdot}
\DeclareFieldFormat[thesis]{title}{\textit{#1\isdot}}
\DeclareFieldFormat[inproceedings]{title}{#1\isdot}
\fi

%-----------------------------------------------------------------------
%maketitle
%-----------------------------------------------------------------------
\newlength{\ifmtitleskip}
\setlength{\ifmtitleskip}{0pt}
\pgfmathsetmacro{\ifm@scaleheader}{12pt/\f@size}
\pgfmathsetlengthmacro{\ifm@headerblockwidth}{0.65\textwidth/\ifm@scaleheader}
\newcommand{\headertitle}{
	%kopfzeile
	\begin{minipage}[b]{40mm}
		\includegraphics[width=33mm]{\@kitlogopath}
	\end{minipage}
	\hfill
	\mbox{\scalebox{\ifm@scaleheader}{%
		\begin{minipage}[b]{\ifm@headerblockwidth}%
				\flushright{\sffamily%
					{\footnotesize%
						\ififmdocument@professorsinheader%
							{\bfseries\normalsize Institut für Mechanik}\\[0.2em]
							Prof.\ Dr.-Ing.\ habil.\ Peter Betsch \\
							Prof.\ Dr.-Ing.\ habil.\ Thomas Seelig
						\else%
							\@headerblock
						\fi%
					}%
				}%
		\end{minipage}
	}}\\
	\vspace{-\baselineskip}
	\vspace{-0.8em}\\
	\rule{\textwidth}{0.2mm}
}
\newcommand{\ifminserttitle}{\Large\sffamily\bfseries\thetitle}
\renewcommand{\maketitle}{%
	\thispagestyle{plain}%
	~\vspace{-\baselineskip}%
	\ififm@usehead\vspace{-\ifm@head}\vspace{-\ifm@headsep}\fi\\
	\headertitle
	\vspace{\ifmtitleskip}
	\begin{center}
		\ifminserttitle
	\end{center}	
	\vspace{0.5em}
}

%----------------------------------------------------------------------
%support for old font commands
%----------------------------------------------------------------------
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}

%----------------------------------------------------------------------
%User commands
%----------------------------------------------------------------------
\newcommand{\aufgabe}[2]{\subsection*{#1. Aufgabe:}{\input{#2}}}%.\dotfill .}}
\newenvironment{mytabular}{\vspace{0.8em}\begin{tabular}{@{}ll}}{\end{tabular}}

\endinput