# 
# d5.g
#
# Evi Haberberger
# August/September 1999
#

sublattice_2_mp := function(Bcon, Bsol, ind, str)

	local i, L, LL, l, ll, z, j, k, s, la, grind, sol, n, lb, lz, a, m;
	
	Print("\nstart sublattice_2_mp():\n\n");
	L := Bsol[1];
	l := Length(L);
	LL := Bcon[2]; 		# conjugacy class of the overgroups
	ll := Length(LL); 
	PrintTo(str, "% file ", str, "\n% created by dode_proto.g\n\n\n");
	AppendTo(str, "input boxes\n\ndefaultfont:=\"cmr7\";\n\n");
	AppendTo(str, "u=0.1mm;\nbeginfig(1);\n\n");
	AppendTo(str, "draw (200u, 0u) -- (200u, 50u);\n");
	AppendTo(str, "draw (200u, 50u) -- (200u, 100u);\n");
	AppendTo(str, "draw (200u, 100u) -- (100u, 300u);\n");
	AppendTo(str, "draw (200u, 100u) -- (300u, 200u);\n");
	AppendTo(str, "draw (300u, 200u) -- (200u, 400u);\n");
	AppendTo(str, "draw (100u, 300u) -- (200u, 400u);\n");
	AppendTo(str, "draw (100u, 300u) -- (0u, 500u);\n");
	AppendTo(str, "draw (200u, 400u) -- (100u, 600u);\n");
	AppendTo(str, "draw (0u, 500u) -- (100u, 600u);\n");
	AppendTo(str, "draw (100u, 600u) -- (200u, 650u);\n");
	AppendTo(str, "circleit.l0(btex 0 etex);\n");
	AppendTo(str, "l0.c = (200u, 0u);\n");
	AppendTo(str, "unfill bpath l0;\n");
	AppendTo(str, "drawboxed(l0);\n");
	AppendTo(str, "circleit.l1(btex 1 etex);\n");
	AppendTo(str, "l1.c = (200u, 50u);\n");
	AppendTo(str, "unfill bpath l1;\n");
	AppendTo(str, "drawboxed(l1);\n");
	AppendTo(str, "circleit.l2(btex 2 etex);\n");
	AppendTo(str, "l2.c = (200u, 100u);\n");
	AppendTo(str, "unfill bpath l2;\n");
	AppendTo(str, "drawboxed(l2);\n");
	AppendTo(str, "circleit.l3(btex 3 etex);\n");
	AppendTo(str, "l3.c = (300u, 200u);\n");
	AppendTo(str, "unfill bpath l3;\n");
	AppendTo(str, "drawboxed(l3);\n");
	AppendTo(str, "circleit.l4(btex 4 etex);\n");
	AppendTo(str, "l4.c = (100u, 300u);\n");
	AppendTo(str, "unfill bpath l4;\n");
	AppendTo(str, "drawboxed(l4);\n");
	AppendTo(str, "circleit.l5(btex 5 etex);\n");
	AppendTo(str, "l5.c = (200u, 400u);\n");
	AppendTo(str, "unfill bpath l5;\n");
	AppendTo(str, "drawboxed(l5);\n");
	AppendTo(str, "circleit.l6(btex 6 etex);\n");
	AppendTo(str, "l6.c = (0u, 500u);\n");
	AppendTo(str, "unfill bpath l6;\n");
	AppendTo(str, "drawboxed(l6);\n");
	AppendTo(str, "circleit.l7(btex 7 etex);\n");
	AppendTo(str, "l7.c = (200u, 650u);\n");
	AppendTo(str, "unfill bpath l7;\n");
	AppendTo(str, "drawboxed(l7);\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[1]));
	AppendTo(str, "$ etex, (180u, 25u));\n");
	AppendTo(str, "label.bot(btex $1$ etex, (180u, 0u));\n");
	AppendTo(str, "label.bot(btex $P$ etex, (180u, 50u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[2]));
	AppendTo(str, "$ etex, (180u, 75u));\n");
	AppendTo(str, "label.bot(btex $N_A(P)$ etex, (180u, 100u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[3]));
	AppendTo(str, "$ etex, (270u, 150u));\n");
	AppendTo(str, "label.bot(btex $A$ etex, (280u, 200u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[4]));
	AppendTo(str, "$ etex, (130u, 200u));\n");
	AppendTo(str, "label.bot(btex $N_{N_G(P)}(A)$ etex, (80u, 300u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[4]));
	AppendTo(str, "$ etex, (230u, 300u));\n");
	AppendTo(str, "label.bot(btex $N_G(A)$ etex, (180u, 400u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[3]));
	AppendTo(str, "$ etex, (170u, 350u));\n");
	AppendTo(str, "label.bot(btex $N_G(P)$ etex, (0u, 520u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(ind[5]));
	AppendTo(str, "$ etex, (30u, 400u));\n");
	AppendTo(str, "label.bot(btex $G$ etex, (180u, 650u));\n");
	AppendTo(str, "label.bot(btex $");
	AppendTo(str, String(Length(Bsol[3][1])));
	AppendTo(str, "$ sol etex, (300u, 190u));\n");
	n := 0;
	a := Size(Bcon[1][1]);		# Bcon[1][1] = A !!
	if l > 1 then
	sol := 0;
	for z in [2..l] do
		j := L[z];
		k := Length(Bcon[2][j]);
		n := n + k;
		s := 250 + (n-k)*20;
		la := 280 + (s - 300)/2;
		m := k mod 2;
		if m = 0 then
			lb := 250 + (n - k/2)*20;
		else
			lb := 250 + (n - (k-1)/2)*20;
		fi; 
		grind := Size(Bcon[1][j])/a;
		AppendTo(str, "draw (");
		AppendTo(str, String(s));
		AppendTo(str, "u, 400u) -- (300u, 200u);\n");
		AppendTo(str, "draw (");
		AppendTo(str, String(s));
		AppendTo(str, "u, 400u) -- (");
		AppendTo(str, String(250 + (n-1)*20));
		AppendTo(str, "u, 400u);\n");
		AppendTo(str, "circleit.l");
		AppendTo(str, String(z+6));
		AppendTo(str, "(btex ");
		AppendTo(str, String(z+6));
		AppendTo(str, " etex);\nl");
		AppendTo(str, String(z+6));
		AppendTo(str, ".c = (");
		AppendTo(str, String(s));
		AppendTo(str, "u, 400u);\n");
		AppendTo(str, "unfill bpath l");
		AppendTo(str, String(z+6));
		AppendTo(str, ";\n");
		AppendTo(str, "drawboxed(l");
		AppendTo(str, String(z+6));
		AppendTo(str, ");\n");
		AppendTo(str, "label.bot(btex $");
		AppendTo(str, String(z));
		AppendTo(str, "$ etex), (");
		AppendTo(str, String(s));
		AppendTo(str, "u, 400u));\n");
		AppendTo(str, "label.bot(btex $");
		AppendTo(str, String(grind));
		AppendTo(str, "$ etex, (");
		AppendTo(str, String(la));
		AppendTo(str, "u, 300u));\n");
		AppendTo(str, "label.bot(btex $");
		AppendTo(str, String(k));
		AppendTo(str, "$ etex, (");
		AppendTo(str, String(lb));
		AppendTo(str, "u, 410u));\n");
		AppendTo(str, "label.bot(btex $");
		AppendTo(str, String(Length(Bsol[3][z])));
		AppendTo(str, "$ sol etex, (");
		AppendTo(str, String(lb));
		AppendTo(str, "u, 390u));\n");
		sol := sol + Length(Bsol[3][z]);
	od;
	fi;
	AppendTo(str, "label.bot(btex in total: ");
	AppendTo(str, String(n+1));
	AppendTo(str, " groups with \\\\");
	AppendTo(str, String(sol));
	AppendTo(str, " solutions etex, (350u, 50u));\n");
	AppendTo(str, "endfig;\n\nend\n");
	Print("\nleave sublattice_2_mp()\n\n");
end;

mp_with_discreta := function(km, strmp)
	
	local cmd, str;
	
	str := "";
	Append(str, String(km));
	Append(str, "_lattice");
	
	cmd := [];
	Append(cmd, "t140e.out ");
	Append(cmd, strmp);	# input filename
	Append(cmd, " ");
	Append(cmd, str);	# output filename
	Append(cmd, " 1");	# manner of output
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

finish_mp_file := function(str)

	local cmd;
	
	cmd := [];
	Append(cmd, "mpost ");
	Append(cmd, str);
	# Append(cmd, ".mp");
	Print(cmd, "\n");
	Exec(cmd);
end;

insert_lattice2tex := function(lat, str, strmp)
	local i;
	
	Print("\nstart insert_lattice2tex():\n\n");
	AppendTo(str, "\n\n\\subsection{Resulting group lattice}\n\n");
	AppendTo(str, "\\input ");
	AppendTo(str, "\\epsfig{figure=");
	AppendTo(str, strmp);
	AppendTo(str, ".1, width=80mm}\n");
	
	# write information concerning the lattice into a tabular:
	# Print("\n\n\n\n", lat, "\n\n\n\n");
	if Length(lat)>2 then
	AppendTo(str, "\\begin{tabular}{|l||c|c|}\n");
	AppendTo(str, "\\hline\nGroup & group order & nb. solutions \\\\ \\hline\n");
		for i in [2..Length(lat)-1] do			# 1=A, Length(lat)=G!
			AppendTo(str, "$B_{");
			AppendTo(str, String(i-1));
			AppendTo(str, "}$ & ");
			AppendTo(str, String(lat[i][1][2]));	# g.o.
			AppendTo(str, " & ");
			AppendTo(str, String(lat[i][1][4]));	# nb. sol.
			AppendTo(str, " \\\\\n");
		od;
	AppendTo(str, "\\hline \n\\end{tabular}\n");
	fi;
	
	Print("\nleave insert_lattice2tex()\n\n");
end;

shadow_latex_report := function(KM, str)
	
	local repstr, texname, l, i, j, strtex;
	
	Print("\nstart shadow_latex_report()\n\n");
	repstr := "report_";
	Append(repstr, str);
	Append(repstr, ".tex");
	PrintTo(repstr, "\\documentclass[11pt]{article}\n");
	AppendTo(repstr, "%titlepage \n");
	AppendTo(repstr, "\\usepackage{fancyheadings} \n");
	AppendTo(repstr, "%\\usepackage{amstex} \n");
	AppendTo(repstr, "\\usepackage{amsmath} \n");
	AppendTo(repstr, "\\usepackage{amssymb} \n");
	AppendTo(repstr, "\\usepackage{multicol} \n");
	AppendTo(repstr, "\\usepackage{epsfig} \n");
	AppendTo(repstr, "\\usepackage{epsf} \n");
	AppendTo(repstr, "\\usepackage{supertabular} \n");
	AppendTo(repstr, "\\usepackage{wrapfig} \n");
	AppendTo(repstr, "%\\usepackage{blackbrd} \n");
	AppendTo(repstr, "%\\usepackage{epic,eepic} \n");
	AppendTo(repstr, "\\usepackage{rotating} \n");
	AppendTo(repstr, "%\\usepackage{concmath} \n");
	AppendTo(repstr, "%\\usepackage{beton} \n\n\n");
	AppendTo(repstr, "%\\pagestyle{empty} \n");
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
	AppendTo(repstr, "\\newcommand{\\EC}{{\\mathbb C}}\n");
	AppendTo(repstr, "\\newcommand{\\cS}{{\\cal S}}\n");
	AppendTo(repstr, "\\newcommand{\\cA}{{\\cal A}}\n");
	AppendTo(repstr, "\\newcommand{\\la}{\\langle}\n");
	AppendTo(repstr, "\\newcommand{\\ra}{\\rangle}\n\n");
	AppendTo(repstr, "\\newcommand{\\eop}{{~~~~\\hspace*{\\fill}{$\\Box$}}}");
	AppendTo(repstr, "\\newcommand{\\eex}{{~~~~\\hspace*{\\fill}{$\\Diamond$}}}");
	AppendTo(repstr, "\\newcommand{\\eod}{{~~~~\\hspace*{\\fill}{$\\diamond$}}}");
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
 	AppendTo(repstr, "   {\\end{trivlist}}\n");
	AppendTo(repstr, "\\newenvironment{proof}%\n");
    	AppendTo(repstr, "    {\\begin{trivlist}%\n");
    	AppendTo(repstr, "    \\item[]{\\it Proof:\\ }}%\n");
    	AppendTo(repstr, "    {\\eop\\end{trivlist}}\n\n\n");
 	AppendTo(repstr, "%\\title{} \n");
	AppendTo(repstr, "%\\author{{\sc Evi Haberberger}} \n");
	AppendTo(repstr, "%\\date{} \n\n");
	AppendTo(repstr, "\\begin{document} \n\n");
	AppendTo(repstr, "%\\maketitle \n\n");
	AppendTo(repstr, "%\\section{Isomorphism classification}\n\n\n");
	Print("\nleave shadow_latex_report()\n\n");
end;

latex_intro := function(strtex)
	Print("\nstart latex_intro():\n\n");
	AppendTo(strtex, "\\section{Theoretical background}\n\n");
	AppendTo(strtex, "The following lemma is easy to verify:\n\n");
	AppendTo(strtex, "\\begin{lemma}\n");
	AppendTo(strtex, "Let $A$ be the (full) stabilizer of two designs\n");
	AppendTo(strtex, "${\\cal D}_1$ and ${\\cal D}_2$, let $g \\in G$ be an\n");
	AppendTo(strtex, "isomorphism between the two designs. Then we have\n");
	AppendTo(strtex, "\\begin{align*}\n");
	AppendTo(strtex, "g \\in N_G(A)\n");
	AppendTo(strtex, "\\end{align*}\n");
	AppendTo(strtex, "\\end{lemma}\n");
	AppendTo(strtex, "\nThis is a necessary condition for the existence of an \n");
	AppendTo(strtex, "isomorphism $g$ between ${\\cal D}_1$ and ${\\cal D}_2$. \n");
	AppendTo(strtex, "So, if we know $A$ to be the \n{\\bf full} automorphism group\n");
	AppendTo(strtex, "of the designs, we only have to look for possible isomorphisms \n");
	AppendTo(strtex, "in the normalizer of $A$ in $G$. But this group quite often\n");
	AppendTo(strtex, "is fairly large and we have to find further restrictions on $g$.\n\n");
	AppendTo(strtex, "On the other hand, it is often not known, if $A$ is the full \n");
	AppendTo(strtex, "automorphism group of the designs.\n");
	AppendTo(strtex, "\\begin{lemma}\\label{NGP}\n");
	AppendTo(strtex, "Let ${\\cal D}_1$ and ${\\cal D}_2$ be two isomorphic designs \n");
	AppendTo(strtex, "on $v$ points and $g \\in G$ an isomorphism between them. \n");
	AppendTo(strtex, "Let $P$ be a $p$-group of $S_v$ such that $P$ is a Sylow \n");
	AppendTo(strtex, "subgroup of $Stab({\\cal D}_1)$ and $Stab({\\cal D}_2)$. \n");
	AppendTo(strtex, "Then there exists a $n \\in N_G(P)$, such that\n\n");
	AppendTo(strtex, "\\begin{align*}\n");
	AppendTo(strtex, "{\\cal D}_1^n = {\\cal D}_2\n");
	AppendTo(strtex, "\\end{align*}\n");
	AppendTo(strtex, "\\end{lemma}\n\n");
	AppendTo(strtex, "Thus, we only need to know at least a suitable $p$-Sylow \n");
	AppendTo(strtex, "subgroup of the stabilizer of a design as a ``good''\n");
	AppendTo(strtex, "guess for the whole stabilizer.\n\n");
	AppendTo(strtex, "\\medskip\n");
	AppendTo(strtex, "When the second lemma is fulfilled, we have two cases:\n");
	AppendTo(strtex, "the first one is $N_G(P) \\leq N_G(A)$; in this case we only \n");
	AppendTo(strtex, "have to form the orbits of $N_G(P)$ on the set of fixed points \n");
	AppendTo(strtex, "of $A$ as $N_G(P)$ acts on this set. In general at least \n");
	AppendTo(strtex, "$N_{N_G(A)}(P)$ acts on the set of fixed points of $A$. \n");
	AppendTo(strtex, "Elements of $N_G(P)$, which are not in $N_{N_G(A)}(P)$ \n");
	AppendTo(strtex, "can be handled by the following lemma:\n\n");
	AppendTo(strtex, "\\begin{lemma}\n");
	AppendTo(strtex, "If there is a $g \\in N_G(P)$ with \n");
	AppendTo(strtex, "$g \\not\\in N_{N_G(A)}(P)$ mapping the design ${\\cal D}_1$ \n");
	AppendTo(strtex, "onto ${\\cal D}_2$, then the group generated by $A$ and its \n");
	AppendTo(strtex, "conjugate $A^g$ is a subgroup of the stabilizer of \n");
	AppendTo(strtex, "${\\cal D}_2$.");
	AppendTo(strtex, "\\end{lemma}\n\n");
	AppendTo(strtex, "After these preparations we now present the theorem \n");
	AppendTo(strtex, "allowing to handle general situations.\n\n");
	AppendTo(strtex, "\\begin{theorem}\n");
	AppendTo(strtex, "Let $G$ act on $\\Omega$, let $A$ a subgroup of $G$ and \n");
	AppendTo(strtex, "$P$ a $p$-Sylow subgroup of $A$.\n");
	AppendTo(strtex, "Let $\\Delta := C_A(\\Omega)$ the set of fixed points of \n");
	AppendTo(strtex, "$A$ and \n\\begin{align*}\n");
	AppendTo(strtex, "\\Gamma = \\bigcup_{g \\in N_G(P), g \\not\\in N_G(A)} \n");
	AppendTo(strtex, "C_{\\la A,A^g \\ra}(\\Omega) \\cup \\bigcup_{A < B \\leq \n");
	AppendTo(strtex, "N_G(A)} C_B(\\Omega)\n");
	AppendTo(strtex, "\\end{align*}\n");
	AppendTo(strtex, "If for $\\delta_1, \\delta_2 \\in \\Delta \\backslash \\Gamma$ \n");
	AppendTo(strtex, "there exists some $g \\in G$ such that $\\delta_1^g = \\delta_2$ \n");
	AppendTo(strtex, "then there exists some $h \\in N_{N_G(P)}(A)$ such that \n");
	AppendTo(strtex, "$\\delta_1^h = \\delta_2$.\n");
	AppendTo(strtex, "\\end{theorem}\n");
	AppendTo(strtex, "\\begin{proof}\n");
	AppendTo(strtex, "We claim first that for each $\\delta \\in \\Delta \\backslash \n");
	AppendTo(strtex, "\\Gamma$ $P$ is a Sylow group of $N_G(\\delta)$.\\\\\n");
	AppendTo(strtex, "Suppose, to the contrary, that for some \n$\\delta \\in \\Delta \n");
	AppendTo(strtex, "\\backslash \\Gamma$ there exists some Sylow group $Q$ with $P < Q$. \n");
	AppendTo(strtex, "Then $P < N_Q(P)$ and there is some $g \\in N_Q(P)$ such that \n");
	AppendTo(strtex, "$P < \\la P,g \\ra$. By the choice of $Q$ and $g \\in Q$, we \n");
	AppendTo(strtex, "have $\\delta^g = \\delta$. So $\\la A,g \\ra$ fixes $\\delta$. \\\\\n");
	AppendTo(strtex, "Because $\\la P,g \\ra > P$, $P$ is Sylow group of $A$ and \n");
	AppendTo(strtex, "$g \\not\\in A$, $A$ is a proper subgroup of $\\la A,g \\ra$. \n");
	AppendTo(strtex, "If $g \\in N_G(A)$ then \n$\\delta \\in \\Gamma$, which is \n");
	AppendTo(strtex, "a contradiction. So $g \\notin N_G(A)$ and $A < \\la A,A^g \n");
	AppendTo(strtex, "\\ra \\leq \\la A,g \\ra$. But then also $\\la A,A^g\\ra$ \n");
	AppendTo(strtex, "fixes $\\delta$, such that again $\\delta \\in \\Gamma$, the \n");
	AppendTo(strtex, "final contradiction. So we can now apply the lemma to complete \n");
	AppendTo(strtex, "the proof.\\end{proof}\n\n");
	AppendTo(strtex, "\\begin{remark}\n");
	AppendTo(strtex, "If $P$ is a $p$-Sylow subgroup of $N_G(A)$, then the proof shows \n");
	AppendTo(strtex, "that \\\\ \n");
	AppendTo(strtex, "$\\Gamma := \\bigcup_{g \\in N_G(P), g \\not\\in N_G(A)} \n");
	AppendTo(strtex, "C_{\\la A,A^g \\ra}(\\Omega)$ suffices in the theorem.\n");
	AppendTo(strtex, "\\end{remark}\n\n");
	AppendTo(strtex, "We apply the following algorithm:\n");
	AppendTo(strtex, "\\begin{enumerate}\n");
	AppendTo(strtex, "\\item Compute $N_G(P)$ and then $N_{N_G(P)}(A)$.\n");
	AppendTo(strtex, "\\item For all $\\delta$ fixed by some overgroup $B \\geq A$ of \n");
	AppendTo(strtex, "$A \\cdot N_{N_G(P)}(A)$ form the same operation and remove these \n");
	AppendTo(strtex, "from $\\Delta$.");
	AppendTo(strtex, "\\item For all $B := \\la A,A^g \\ra$, $g$ elements of a transversal \n ");
	AppendTo(strtex, "of $N_G(P)$ over $N_{N_G(P)}(A)$, remove the fixed points of $B$ \n");
	AppendTo(strtex, "from $\\Delta$; handle these fixed points for $B$ similarly.\n");
	AppendTo(strtex, "\\item Determine the orbits of $N_{N_G(P)}(A)$ on the remaining \n");
	AppendTo(strtex, "set $\\Delta$.\n");
	AppendTo(strtex, "\\end{enumerate}\n");
	AppendTo(strtex, "The construction of $\\langle A, A^g \\rangle$ quite often gives\n");
	AppendTo(strtex, "the symmetric or alternating group on $v$ points. These groups \n");
	AppendTo(strtex, "cannot arise as automorphism groups of nontrivial designs, so \n");
	AppendTo(strtex, "we can exclude these cases. Another fact is, that we always get a\n");
	AppendTo(strtex, "better approach to the full automorphism group of the designs\n");
	AppendTo(strtex, "with this construction.\n\n");
	AppendTo(strtex, "%\\medskip\n");
	AppendTo(strtex, "%In order to be able to apply lemma \\ref{NGP} we only\n");
	AppendTo(strtex, "%consider the groups with $P$ as Sylow subgroup.\n");
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

