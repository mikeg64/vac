#!/usr/bin/perl

# Usage: html2latex HTMLFILE > TEXFILE

# Use this package to replace TAB-s with spaces
use Text::Tabs;

$tabstop=8;

# Horizontal line in verbatim environment
$hr= '-' x 79 . "\n";

while(<>){

    # Special HTML characters: &amp; <-- & , &gt; <-- > , &lt; <-- <
    s/\&/\&amp;/g;
    s/\$>\$/\&gt;/g;
    s/\$<\$/\&lt;/g;
    s/>/\&gt;/g;
    s/</\&lt;/g;

    # <PRE> <-- \begin{verbatim}
    $verbatim=1 if s/\\begin\{verbatim\}/<pre>/ig;

    if($verbatim){
	# </PRE> <-- \end{verbatim}
	$verbatim=0 if s/\\end\{verbatim\}/<\/pre>/ig;

        # Use <hr> for horizontal line in verbatim environment
	s/$hr/<hr>/i;

	# Replace TABS by spaces and print line
	print expand($_);

	next;
    }

    # Paragraph added for two empty lines
    if(/^$/){
	if($prevempty){
	    $prevempty = 0;
	    print "<p>\n";
	}else{
	    $prevempty = 1;
	}
	next;
    }
    $prev_empty = 0;

    # Replace special Tex characters

    # \ and $: \ <-- \backslash  ; $ <-- \$ ; \backslash <-- $\backslash$
    s/\$\\backslash\$/\\backslash/g;
    s/\\\$/\$/g;
    s/\\backslash/\\/g;

    # Brackets: [ <-- $[$ , ] <-- $]$ , { <-- $\{$ , } <-- $\}$ , | <-- $|$
    s/\$\[\$/\[/g;
    s/\$\]\$/\]/g;
    s/\$\\\{\$/\{/g;
    s/\$\\\}\$/\}/g;
    s/\$\|\$/\|/g;
    s/\$\+\$/<b>\+<\/b>/g;
    s/\$\-\$/<b>\-<\/b>/g;
    s/\$([\-\+\*])\$/$1/g;

    # ^,~,_,%: ^ <-- \^{} , ~ <-- \~{} , _ <-- \_ , % <-- \% , # <-- \#
    s/\\^{}/\^/g;
    s/\\~{}/\~/g;
    s/\\_/\_/g;
    s/\\%/\%/g;
    s/\\#/\#/g;

    # Translate HTML accents to Tex: 
    # &.acute; --> \'. , &.uml; --> \". , &.circ; --> \~.
    s/\\'(.)/\&($1)acute;/g;
    s/\\"(.)/\&($1)uml;/g;
    s/\\~(.)/\&($1)circ;/g;

    # Replace "This section:" with "This page:"
    s/^This section:/This page:/;

    # Replace \\newpage with </HTML>
    s/\\newpage/<\/HTML>/;

    # Move IMG out of title
    # s/(<h\d>)(<IMG [^>]*> *)/$2$1/i;

    # Single formatting symbols
    # s/^$/<p>/g;
    s/\\\\/<br>/g;
    s/\\medskip\\hrule\\medskip /<hr>/g;

    # Beginning of HTML font styles <-- {\TexFont
    s/\{\\bf ([^}]+)\}/<b>$1<\/b>/ig;
    s/\{\\it ([^}]+)\}/<i>$1<\/i>/ig;
    s/\{\\em ([^}]+)\}/<em>$1<\/em>/ig;
    s/\{\\tt ([^}]+)\}/<kbd>$1<\/kbd>/ig;

    # Beginning of titles --> \TEXSection{
    s/\\section\{([^\}]+)\}/<H1>$1<\/H1>/g;
    s/\\subsection\{([^\}]+)\}/<H2>$1<\/H2>/g;
    s/\\subsubsection\{([^\}]+)\}/<H4>$1<\/H4>/g;

    # Centering
    s/\\begin{center}/<center>/ig;
    s/\\end{center}/<\/center>/ig;

    # Lists
    s/\\begin{itemize}/<ul>/ig;
    s/\\item /<li>/ig;
    s/\\end{itemize}/<\/ul>/ig;

    s/\\begin{enumerate}/<ol>/ig;
    s/\\end{enumerate}/<\/ol>/ig;

    s/<dl>/\\begin{description}/ig;
    s/\\item\[([^\]]*)\]/<dt>$1<\/dt><dd>/ig;
    s/\\end{description}/<\/dl>/ig;

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

    # Replace the remaining special characters

    #s/</\$<\$/g;
    #s/>/\$>\$/g;
    #s/\&/\\\&/g;

    # Print line
    print;
}

exit;
