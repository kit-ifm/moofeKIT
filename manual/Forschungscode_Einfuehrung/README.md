# Vorlage für allgemeine Dokumente am IFM
Soll für allgemeine Dokumente verwendet werden. Das Hauptpaket [ifmdocument](ifmdocument.sty) baut auf scrartcl auf und ist Basis für einige weitere Vorlagen:
* [BA_MA_Aufgabenstellung](./BA_MA_Aufgabenstellung/)
* [BA_MA_Bewertung](./BA_MA_Bewertung/)
* [Uebungsblatt] (./Uebungsblatt/)

Die **Hauptdatei** ist [main.tex](main.tex).
Dort werden wichtige Basisinfos angegeben und die Dokumentstruktur
(aus mehreren files) aufgebaut. Im [header](header.tex) sind bereits enige häufig
gebrauchte Befehle definiert. Diese können nach belieben ergänzt werden. 
Zusätzliche Pakete und Befehlsdefinitionen sollten im [header](header.tex)
eingefügt werden. 

## Optionen des Package ifmdocument
Alle Optionen der KOMAScript-Klasse scrartcl sind verfügbar und werden an diese Klasse weitergereicht. 
### Allgemeines
- **professorsinheader:** falls true werden die professoren im header der ersten
	Seite rechts gedruckt. Sollte nur bei offiziellen, von beiden Professoren 
	abgesegneten, Dokumenten ausgewählt werden. Normalfall: Befehl \headerblock{} nutzen. 
- **language:** 'english' oder 'german' (standard). Dokumentsprache. z.b. wird automatisch das richtige KIT-Logo ausgewählt. 
- **10pt,12pt,...:** Schriftgröße. Empfohlen 12pt. 
- **twoside/oneside:** ein- oder zweiseitiger Druck (Duplex).
### Kopf- und Fußzeile
Die Kopfzeile wird auf der ersten Seite nicht gesetzt. Dort wird eine Kopfzeile mit dem KIT-Logo und informationen aus dem headerblock bzw. professorsinheader gedruckt.
- **pagenumberinfooter:** setzt die Seitenzahl in der Fußzeile (außen)
- **pagenumberinheader:** setzt die Seizenzahl in der Kopfzeile (außen) und aktivierte den header. 
- **sectioninheader:** setzt das aktuelle Kapitel in die Kopfzeile (innen)
- **authorinfooter:** setzt den mit \author{} definierten Author in die Fußzeile (innen)
- **dateinfooter:** setzt das mit \date{} definierte Datum in die Fußzeile (innen), kann mit 
	"authorinfooter" kombiniert werden 
### Bibliographie
- **bibbackend:** definiert welches Bibliographieprogramm (biber oder bibtex)
 	verwendet wird. Empfohlen wird **bibtex**. Dies muss auch in den
	Einstellungen des TeX-Programms gewählt werden (Kompilierungskette)
- **bibstyle:** Standard ist 'none'. Verfügbare Zitationsstile alpha \[ABC92\], authoryear oder numeric [1]. Empfohlen
	ist **numeric** aber insbesondere während der Erstellung bietet alpha Vorteile, da es einfacher lesbar ist.

## Grundlegende Angaben
Angaben wie Author, titel etc. werden in der Präembel des Latex-Dokumentes gesetzt
(häufig im [header](header.tex). Besondere Befehle für die IFM-Klasse sind:
- **\kitlogo{pfad}:** pfad zum KITlogo. Standard ist "kitlogo" bzw. "kitlogoEN" für englisch-sparchige Texte (language=english). 
- **\headerblock{header}:** Inhalt des headers rechts auf der ersten Seite. Kann 
	mit beliebigen Formatierungsoptionen gefüllt werden. 

## Userbefehle, -pakete
Neue Befehle und zusätzliche Pakete sollen in [header.tex](header.tex)
eingefügt werden.

## Befehlsvervollständigung
Damit die Befehle aus [ifmdocument](ifmdocument.cls) richtig erkannt werden,
müssen diese häufig dem Texteditor bekannt gemacht werden.
- **Texstudio**
	- Windows+MikTex: kopiere die datei [ifmdocument.cwl](ifmdocument.cwl) in
	den Ordner
	C:\Users\DeinUsername\AppData\Roaming\texstudio\completion\user
	
## Tikz Externalization
Es besteht die Möglichkeit, tikz-Bilder als eigenständige files (eps/pdf) erzeugen zu lassen, damit das dokument schneller kompiliert werden kann. 
Dazu muss "--shell-escape" als option beim Aufrufen des latex/pdflatex befehls übergeben werden. 
- **TeXstudio Var1:** Magic-Comment in der Hauptdatei <br>
*%!TeX TXS-program:latex = latex.exe -src -interaction=nonstopmode --shell-escape %.tex* bzw. <br>
*%!TeX TXS-program:pdflatex = pdflatex -synctex=1 -interaction=nonstopmode --shell-escape %.tex*  
- **TeXstudio Var2:** Eintragen der Option "--shell-escape" unter Optionen->TeXStudio konfigurieren -> Befehle -> LaTeX/PdfLaTeX (Achtung: Variante 2 stellt ein Sicherheitsrisiko dar, da potentiell jedes Dokument, das kompiliert wird, vollen Systemzugriff erhält.)
- **TeXmaker:** Die Varianten unter TeXstudio sollten hier auch funktionieren (nicht getestet)