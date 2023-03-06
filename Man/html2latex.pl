#!/usr/bin/perl

# Usage: html2latex HTMLFILE > TEXFILE

# Use this package to replace TAB-s with spaces
use Text::Tabs;

$tabstop=8;

# Horizontal line in verbatim environment
$hr= '-' x 79 . "\n";

# The beginning of the LATEX document
print "\\input MANhead.tex\n";


while(<>){

    # Counter for <!-- LATEX commands --> in one line
    $ilatex=0;

    # Hide the LATEX commands into ' MYLATEX0 ', ' MYLATEX1 ', ' MYLATEX2 ' ...
    # and save content into @latex
    while(s/<\!\-\- LATEX (.*?) \-\->/ MYLATEX$ilatex /i){
	$latex[$ilatex++]=$1;
    }
 
    # Comment out HTML header and ADDRESS
    $comment=1 if(s/<(HEAD|ADDRESS)>//);
    if($comment){	
	$comment=0 if(s/<\/(HEAD|ADDRESS)>//);
	print "\%$_";
	next;
    }

    # <PRE> --> \begin{verbatim}
    $verbatim=1 if s/<pre>/\\begin\{verbatim\}/ig;

    if($verbatim){
	# </PRE> --> \end{verbatim}
	$verbatim=0 if s/<\/pre>/\\end\{verbatim\}/ig;

        # Use $hr for horizontal line in verbatim environment
	s/<hr[^>]*>/$hr/i;

	# Remove all other HTML tags
	s/<[^>]+>//gi;

        # Get rid off accents in verbatim environment
	s/\&(.)(acute|uml|circ);/$1/g;

	# Special HTML characters: &amp; --> & , &gt; --> > , &lt; --> <
	s/\&amp;/\&/g;
	s/\&gt;/>/g;
	s/\&lt;/</g;

	# Replace TABS by spaces and print line
	print expand($_);

	next;
    }

    # Lines broken in the middle of <...> are put together
    if(/<[^>]*$/){
	chop;
	$cont=$_;
	next;
    }
    if($cont){
	$_=$cont.$_;
	$cont=0;
    }

    # Replace special Tex characters

    # \ and $: \ --> \backslash  ; $ --> \$ ; \backslash --> $\backslash$
    s/\\/\\backslash/g;
    s/\$/\\\$/g;
    s/\\backslash/\$\\backslash\$/g;

    # Brackets: [ --> $[$ , ] --> $]$ , { --> $\{$ , } --> $\}$ , | --> $|$
    s/\[/\$[\$/g;
    s/\]/\$]\$/g;
    s/\{/\$\\\{\$/g;
    s/\}/\$\\\}\$/g;
    s/\|/\$\|\$/g;
    s/<b>\+<\/b>/\$\+\$/g;
    s/<b>\-<\/b>/\$\-\$/g;
    s/([\-\+\*])(\d)/\$$1\$$2/g;
    s/(\d)([\-\+\*])/$1\$$2\$/g;

    # $$ --> ''  e.g. {} --> ${$$}$ --> ${}$
    s/(?!\\)\$\$//g;
    
    # ^,~,_,%: ^ --> \^{} , ~ --> \~{} , _ --> \_ , % --> \% , # --> \#
    s/\^/\\^{}/g;
    s/\~/\\~{}/g;
    s/\_/\\_/g;
    s/\%/\\%/g;
    s/\#/\\#/g;

    
    # Translate HTML accents to Tex: 
    # &.acute; --> \'. , &.uml; --> \". , &.circ; --> \~.
    s/\&(.)acute;/\\'$1/g;
    s/\&(.)uml;/\\"$1/g;
    s/\&(.)circ;/\\~$1/g;

    # Replace anchors
    while(/<(A |\/A>)/){
       $http=$1 if s/<A HREF=\"?(http:\/\/[^>"]*)\"?>/{\\bf /i;
       $href=1 if s/<A HREF=[^>]*>/{\\bf /i;
       s/<A[^>]*>//gi;
       if($http){
           $http=0 if s/<\/A>/} \($http\)/i;
       }elsif($href){
           $href=0 if s/<\/A>/}/i;
       }else{
           s/<\/A>//i;
       }
    }

    # Remove some lines completely
    # We don't want the vaclogo repeated inside the manual
    s/<IMG SRC="vaclogo.gif"[^>]*><br>//i;
    s/<i>Logo by G. Egedi<\/i>//;
 
    s/<\/?BODY[^>]*>//;
    s/<HTML>//;

    # Replace "This page:" with "This section:"
    s/^This page:/This section:/;

    # Start new page at the end of each HTML file
    s/<\/HTML>/\\newpage/;

    # Move IMG out of title
    s/(<h\d>)(<IMG [^>]*> *)/$2$1/i;

    # Remove ordered lists from titles
    $dropol4=1 if s/<ol><h4/<h4/i;
    $dropol3=1 if s/<ol><h3/<h3/i;
    s/<li>//i if /<h\d>/i;
    $dropol4=0 if $dropol4 && s/<\/ol>//i;
    $dropol3=0 if $dropol3 && s/<\/ol>//i;

    # Special HTML characters: &amp; --> & , &gt; --> > , &lt; --> <
    s/\&amp;/\&/g;
    s/\&gt;/>/g;
    s/\&lt;/</g;

    # Single formatting symbols
    s/<p>/\n/g;
    s/<br>/\\\\/ig;
    s/<hr[^>]*>/\\medskip\\hrule\\medskip /ig;

    # Beginning of HTML font styles  --> {\TexFont
    s/<b>/\{\\bf /ig;
    s/<i>/\{\\it /ig;
    s/<em>/\{\\em /ig;
    s/<kbd>/\{\\tt /ig;
    s/<strong>/\{\\bf /ig;

    # Beginning of titles --> \TEXSection{
    s/<h1> */\\section\{/ig;
    s/<h2> */\\subsection\{/ig;
    s/<h3> */\\subsection\{/ig;
    s/<h4> */\\subsubsection\{/ig;

    # End of font styles and titles --> }
    s/<\/(b|i|em|kbd|strong|h1|h2|h3|h4)>/\}/ig;

    # Centering
    s/<center>/\\begin{center}/ig;
    s/<\/center>/\\end{center}/ig;

    # Lists
    s/<ul>/\\begin{itemize}/ig;
    s/<li>/\\item /ig;
    s/<\/ul>/\\end{itemize}/ig;

    s/<ol>/\\begin{enumerate}/ig;
    s/<li>/\\item /ig;
    s/<\/ol>/\\end{enumerate}/ig;

    s/<dl>/\\begin{description}/ig;
    s/<dt>/\\item[/ig;
    s/<\/dt>//ig;
    s/<dd>/]/ig;
    s/<\/dl>/\\end{description}/ig;

    # GIF Images 
    while(/<IMG SRC=\"?([\w\.]*?)\.gif/){
	$name=$1;
	if(-f "$name.tex"){
            # <IMG SRC=filename.gif> --> \input filename.tex
	    s/<IMG SRC=[^>]*>/\\input $name.tex\n/;
	}else{
            # <IMG SRC=filename.gif> --> \epsfig{file=filename.eps}
	    s/<IMG SRC=[^>]*>/\\epsfig{file=$name.eps}/;

	    # Use convert to do the conversion if necessary
	    `convert $name.gif $name.eps` unless -e "$name.eps";
	}
    }

    # HTML FORM
    s/<FORM>/\n\\bigskip\\noindent/;
    s/<\/FORM>/\n\\bigskip\n/;
    s/<INPUT TYPE=checkbox>/ \\framebox[1em]{\\vphantom{x}} /gi;
    s/<INPUT TYPE="?submit"? VALUE=(\w+)[^>]*>/ \\fbox{$1} /gi;
    s/<INPUT TYPE="?submit"? VALUE="([^"]+)[^>]*>/ \\fbox{$1} /gi;

    if(s/<INPUT TYPE=radio([^>]*)>//i){
        $head=$`; $tail=$';
        $x = ($1=~/CHECKED/)? "x" : '\\vphantom{x}';
        $_=$head." \\framebox[1em]{$x} ".$tail;
    }

    if(s/<INPUT ([^>]*)>//i){
        $head=$`; $tail=$'; $arg=$1;
        $size = ($arg=~/SIZE=(\d+)/)? $1 : 15;
        $_=$head." \\framebox[$size em]{\\vphantom{T}} ".$tail;
    }

    s/<SELECT><OPTION SELECTED>([^<]*).*?<\/SELECT>/ \\fbox{$1} /i;

    # Replace the remaining special characters

    s/</\$<\$/g;
    s/>/\$>\$/g;
    s/\&/\\\&/g;

    # Insert back the Latex text: ' MYLATEX12 ' --> $latex[12]

    s/ MYLATEX(\d+) /$latex[$1]/g;

    # Print line
    print;
}

# Finish LATEX document
print "\\end{document}\n";

exit;
