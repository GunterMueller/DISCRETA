# i3.g
#
# Evi Haberberger
# December 1999
#


finish_mp_file := function(str)

	local cmd;
	
	cmd := [];
	Append(cmd, "mpost ");
	Append(cmd, str);
	# Append(cmd, ".mp");
	Print(cmd, "\n");
	Exec(cmd);
end;

mp_with_discreta := function(km, strmp)
	
	local cmd, str;
	
	str := "lattice_";
	Append(str, String(km));
	
	cmd := [];
	Append(cmd, "t140e.out mpbase_");
	Append(cmd, String(km));
	Append(cmd, " ");
	Append(cmd, str);
	Append(cmd, " 1");
	Print(cmd, "\n");
	Exec(cmd);
	
	# mpost the resulting file
	cmd := [];
	Append(cmd, "mpost ");
	Append(cmd, str);
	Append(cmd, ".mp");
	Print(cmd, "\n");
	Exec(cmd);
end;

insert_lattice2tex := function(str, strmp)
	Print("\nstart insert_lattice2tex():\n\n");
	AppendTo(str, "\n\n\\subsection{Resulting group lattice}\n\n");
	AppendTo(str, "\\input ");
	AppendTo(str, "\\epsfig{figure=./");
	AppendTo(str, strmp);
	AppendTo(str, ".1, width=120mm}\n");
	Print("\nleave insert_lattice2tex()\n\n");
end;

shadow_latex_report := function(KM, str)
	
	local repstr, texname, l, i, j, strtex;
	
	Print("\nstart shadow_latex_report()\n\n");
	repstr := "report_";
	Append(repstr, str);
	Append(repstr, ".tex");
	PrintTo(repstr, "\\documentclass[]{article}\n");
	AppendTo(repstr, "%12pt,titlepage \n");
	AppendTo(repstr, "\\usepackage{fancyheadings} \n");
	AppendTo(repstr, "%\\usepackage{amstex} \n");
	AppendTo(repstr, "\\usepackage{amsmath} \n");
	AppendTo(repstr, "\\usepackage{amssymb} \n");
	AppendTo(repstr, "\\usepackage{multicol} \n");
	AppendTo(repstr, "\\usepackage{epsfig} \n");
	AppendTo(repstr, "\\usepackage{supertabular} \n");
	AppendTo(repstr, "\\usepackage{wrapfig} \n");
	AppendTo(repstr, "%\\usepackage{blackbrd} \n");
	AppendTo(repstr, "%\\usepackage{epic,eepic} \n");
	AppendTo(repstr, "\\usepackage{rotating} \n");
	AppendTo(repstr, "%\\usepackage{concmath} \n");
	AppendTo(repstr, "%\\usepackage{beton} \n\n\n");
	AppendTo(repstr, "\\pagestyle{empty} \n");
	AppendTo(repstr, "\\evensidemargin 0in \n");
	AppendTo(repstr, "\\oddsidemargin 0in \n");
	AppendTo(repstr, "\\marginparwidth 0pt \n");
	AppendTo(repstr, "\\marginparsep 0pt \n\n");
	AppendTo(repstr, "\\topmargin -1in \n");
	AppendTo(repstr, "\\headheight 0.7cm \n");
	AppendTo(repstr, "\\headsep 1.8cm \n");
	AppendTo(repstr, "%\\footheight 0.7cm \n");
	AppendTo(repstr, "\\footskip 2cm \n");
	AppendTo(repstr, "\\textheight 22cm \n");
	AppendTo(repstr, "\\textwidth 6.2in \n\n\n");
	AppendTo(repstr, "\\marginparpush 0pt \n\n\n\n");
	AppendTo(repstr, "\\renewcommand{\\baselinestretch}{1.5} \n\n");
	AppendTo(repstr, "\\newcommand{\\EN}{{\\mathbb N}}\n");
	AppendTo(repstr, "\\newcommand{\\EQ}{{\\mathbb Q}}\n");
	AppendTo(repstr, "\\newcommand{\\ER}{{\\mathbb R}}\n");
	AppendTo(repstr, "\\newcommand{\\EZ}{{\\mathbb Z}}\n");
	AppendTo(repstr, "\\newcommand{\\EC}{{\\mathbb C}}\n\n");
	AppendTo(repstr, "\\newcounter{theoremcounter}[section]\n");
	AppendTo(repstr, "\\def\\thetheoremcounter{\\thesection.\\arabic{theoremcounter}}\n");
	AppendTo(repstr, "\\newcommand{\\labell}[1]{\\label{#1}%\n");
	AppendTo(repstr, "\\ifmmode $$\\vspace*{-\\baselineskip}\\marginpar{#1}%\n");
	AppendTo(repstr, "\\vspace*{-\\baselineskip}$$\\else\\marginpar{#1}\\fi}\n");
	AppendTo(repstr, "\\newenvironment{theorem}%\n");
	AppendTo(repstr, "    {\\begin{trivlist}\\refstepcounter{theoremcounter}%\n");
	AppendTo(repstr, "    \\item[]{\\bf\\thetheoremcounter\\ Theorem\\ }\\em}%\n");
	AppendTo(repstr, "    {\\end{trivlist}}\n");
	AppendTo(repstr, "\\newenvironment{definition}%\n");
	AppendTo(repstr, "    {\\begin{trivlist}\\refstepcounter{theoremcounter}%\n");
	AppendTo(repstr, "    \\item[]{\\bf\\thetheoremcounter\\ Definition\\ }}%\n");
	AppendTo(repstr, "    {\\end{trivlist}}\n");
	AppendTo(repstr, "\\newenvironment{corollary}%\n");
	AppendTo(repstr, "    {\\begin{trivlist}\\refstepcounter{theoremcounter}%\n");
	AppendTo(repstr, "    \\item[]{\\bf\\thetheoremcounter\\ Corollary\\ }\\em}%\n");
	AppendTo(repstr, "    {\\end{trivlist}}\n");
	AppendTo(repstr, "\\newenvironment{defundsatz}%\n");
	AppendTo(repstr, "    {\\begin{trivlist}\\refstepcounter{theoremcounter}%\n");
	AppendTo(repstr, "    \\item[]{\\bf\\thetheoremcounter\\ Definition and theorem\\ }\\em}%\n");
	AppendTo(repstr, "    {\\end{trivlist}}\n");
	AppendTo(repstr, "\\newenvironment{remark}%\n");
	AppendTo(repstr, "    {\\begin{trivlist}\\refstepcounter{theoremcounter}%\n");
	AppendTo(repstr, "    \\item[]{\\bf\\thetheoremcounter\\ Remark\\ }}%\n");
	AppendTo(repstr, "    {\\end{trivlist}}\n");
	AppendTo(repstr, "\\newenvironment{lemma}%\n");
	AppendTo(repstr, "    {\\begin{trivlist}\\refstepcounter{theoremcounter}%\n");
	AppendTo(repstr, "    \\item[]{\\bf\\thetheoremcounter\\ Lemma\\ }\\em}%\n");
 	AppendTo(repstr, "   {\\end{trivlist}}\n\n\n");
	AppendTo(repstr, "%\\title{} \n");
	AppendTo(repstr, "%\\author{{\sc Evi Haberberger}} \n");
	AppendTo(repstr, "%\\date{} \n\n");
	AppendTo(repstr, "\\begin{document} \n\n");
	AppendTo(repstr, "%\\maketitle \n\n");
	AppendTo(repstr, "\\section{Isomorphism classification}\n\n\n");
	Print("\nleave shadow_latex_report()\n\n");
end;

latex_intro := function(strtex)
	Print("\nstart latex_intro():\n\n");
	AppendTo(strtex, "\\subsection{Theoretical background}\n\n");
	AppendTo(strtex, "If we want to find isomorphisms between $t-(v, k, \\lambda)$ designs \n");
	AppendTo(strtex, "we already constructed, there is a very important\n");
	AppendTo(strtex, "and basic lemma we can use:\n\n");
	AppendTo(strtex, "\\begin{lemma}\n");
	AppendTo(strtex, "Let $A$ be the (full) stabilizer of two designs\n");
	AppendTo(strtex, "${\\cal D}_1$ and ${\\cal D}_2$, let $g \\in G$ be an\n");
	AppendTo(strtex, "isomorphism between the two designs. Then we have\n");
	AppendTo(strtex, "\\begin{align*}\n");
	AppendTo(strtex, "g \\in N_G(A)\n");
	AppendTo(strtex, "\\end{align*}\n");
	AppendTo(strtex, "\\end{lemma}\n");
	AppendTo(strtex, "So, if we know $A$ to be the full automorphism group\n");
	AppendTo(strtex, "of a design, we only have to look for possible isomorphisms \n");
	AppendTo(strtex, "in the normalizer of $A$ in $G$. But this group quite often\n");
	AppendTo(strtex, "is fairly large and we have to find further restrictions to $g$.\n\n");
	AppendTo(strtex, "\\begin{lemma}\\label{NGP}\n");
	AppendTo(strtex, "Let ${\\cal D}_1$ and ${\\cal D}_2$ be two isomorphic designs on $v$ points\n");
	AppendTo(strtex, "and $g \\in G$ the isomorphism between them. Let $P$ be a $p$-group of $S_v$\n");
	AppendTo(strtex, "such that $P$ is subgroup of $Stab({\\cal D}_1)$ and Sylow subgroup of \n");
	AppendTo(strtex, "$Stab({\\cal D}_2)$. Then we can find a $n \\in N_G(P)$ already, such that\n\n");
	AppendTo(strtex, "\\begin{align*}\n");
	AppendTo(strtex, "{\\cal D}_1^n = {\\cal D}_2\n");
	AppendTo(strtex, "\\end{align*}\n");
	AppendTo(strtex, "\\end{lemma}\n\n");
	AppendTo(strtex, "Thus, we need to know a suitable $p$-Sylow subgroup\n");
	AppendTo(strtex, "of the stabilizer of a design - or at least of a ``good''\n");
	AppendTo(strtex, "guess for this stabilizer.\n\n");
	AppendTo(strtex, "\\medskip\n");
	AppendTo(strtex, "In case the second lemma is fulfilled, we have two cases:\n");
	AppendTo(strtex, "the first one is $A \\leq N_G(P)$; in this case we only have to\n");
	AppendTo(strtex, "form the orbits of $N_G(P)$ on the set of fixed points of $A$\n");
	AppendTo(strtex, "as $N_G(P)$ acts on this set. In case $A$ is no subgroup of\n");
	AppendTo(strtex, "$N_G(P)$ it is easy to see that $N_{N_G(A)}(P)$ acts on \n");
	AppendTo(strtex, "the set of fixed points of $A$. Then we have to use \n");
	AppendTo(strtex, "another method according to the following lemma:\n\n");
	AppendTo(strtex, "\\begin{lemma}\n");
	AppendTo(strtex, "If there is a $g \\in N_G(P)$ with $g \\not\\in N_{N_G(A)}(P)$\n");
	AppendTo(strtex, "mapping the design ${\\cal D}_1$ isomorphically on ${\\cal D}_2$,\n");
	AppendTo(strtex, "then the group generated by $A$ and its conjugate $A^g$ is a subgroup\n");
	AppendTo(strtex, "of the stabilizer of ${\\cal D}_2$.\n");
	AppendTo(strtex, "\\end{lemma}\n\n");
	AppendTo(strtex, "The construction of $\\langle A, A^g \\rangle$ quite often gives\n");
	AppendTo(strtex, "the symmetric or alternating group on $v$ points. These groups cannot arise as\n");
	AppendTo(strtex, "automorphism groups of nontrivial designs, so we can exclude\n");
	AppendTo(strtex, "these cases. A very interesting fact is, that we always get a\n");
	AppendTo(strtex, "better approach to the full automorphism group of the designs\n");
	AppendTo(strtex, "with this construction.\n\n");
	AppendTo(strtex, "\\medskip\n");
	AppendTo(strtex, "In order to be able to apply lemma \\ref{NGP} we only\n");
	AppendTo(strtex, "consider those groups having $P$ as Sylow subgroup.\n");
	AppendTo(strtex, "\n");
	AppendTo(strtex, "\n");
	AppendTo(strtex, "\n");
	Print("\nleave latex_intro()\n\n");
end;

latex_results := function(str)
	
	local cmd, str1;
	
	Print("\nstart latex_results():\n\n");
	str1 := [];
	Append(str1, str);
	Append(str1, ".tex");
	AppendTo(str1, "\n\n\\end{document} \n");
	
	cmd := [];
	Append(cmd, "latex ");
	Append(cmd, str1);
	Print(cmd, "\n");
	Exec(cmd);
	Print(cmd, "\n");
	Exec(cmd);
	Print("\nleave latex_results()\n\n");
end;

create_ps_file_and_view := function(str)

	local cmd, gh;
	
	cmd := [];
	Append(cmd, "dvips ");
	Append(cmd, str);
	Append(cmd, ".dvi -o");
	Print(cmd, "\n");
	Exec(cmd);
	
	gh := [];
	Append(gh, "ghostview ");
	Append(gh, str);
	Append(gh, ".ps &");
	Print(gh, "\n");
	Exec(gh);
end;

