# d3.g
#
# Evi Haberberger
# August/September 1999
#


compute_fuse := function(solcon, KM, Bsol, Bcon, param, strtex, fv) 
# Now the fusing process can be started:
# compute fusing of orbits of A under the relevant groups in Borb[]:
# Fuse[] is fusion mapping 

	local Fuse, len, j, len1, fuse, i, p, str1, f, k, l, n, z, ind, ll, x;

	Print("\n\nstart compute_fuse():\n");
	if fv > 0 then
		AppendTo(strtex, "\\subsection{Fusion process and ``true solutions''}\n\n");
		AppendTo(strtex, "Compute the fusion mapping of $A$ under all groups with solutions:\n");
		if fv > 2 then
			if Length(Bsol[1]) > 1 then
				AppendTo(strtex, "In the following matrix the rows correspond to the orbits of $A$ on the\n");
				AppendTo(strtex, String(param[3]));
				AppendTo(strtex, "-subsets of our point set. The columns correspond to the above constructed\n");
				AppendTo(strtex, "overgroups (including $A$ in the first column). The entry $M_{i, j}$\n");
				AppendTo(strtex, "is the index of the image orbit of orbit $i$ (w.r.t. A) under the action of\n");
				AppendTo(strtex, "group $j$. The first row gives the identity mapping, because\n");
				AppendTo(strtex, "the orbits of $A$ already were canonical.\n");
			fi;
		fi;
		AppendTo(strtex, "\n\\smallskip\n");
	fi;
	Fuse := [];
	len := Length(solcon[1][1][1]);
	for j in [1..Length(Bsol[1])] do
		len1 := Length(solcon[2][j]);
		fuse := [];
		for i in [1..len1] do
#			p := PermList(solcon[2][j][i]);
			str1 := StructuralCopy(KM);
			Append(str1, "_");
			Append(str1, String(Bsol[1][j]));
			Append(str1, "_");
			Append(str1, String(i));
			f := fuse_orbits_by_representatives(solcon[1][1][1], str1);
		  	Add(fuse, f[2]);
		od;
		Add(Fuse, fuse);
	od;
	
	if fv > 2 and Length(solcon[1][1][1]) < 200 then
		if Length(Bsol[1]) > 1 then
			# print to file strtex:
			l := Length(solcon[1]);
			n := Length(solcon[1][1][1]);
			AppendTo(strtex, "{\\tiny\n\\begin{supertabular}{|");
			for z in [1..l] do
				ind := Bsol[1][z];
				AppendTo(strtex, "*{");
				AppendTo(strtex, String(Length(Bcon[2][ind])));
				AppendTo(strtex, "}{c}");
				AppendTo(strtex, "|");
			od;
			AppendTo(strtex, "}\\hline\n");
			for i in [1..n] do
				x := 0;
				for z in [1..l] do
					ll := Length(solcon[1][z]);
				#	Print("x=", x, " ll=", ll, "\n");
				#	Print("Fuse=", Fuse, "\n");
					for j in [1..ll] do
				#		Print("x+j=", x+j, " i=", i, "\n");
						if Fuse[z][j][i] <> 0 then
							AppendTo(strtex, "{\\bf ");
						fi;
						AppendTo(strtex, Fuse[z][j][i]);
						if Fuse[z][j][i] <> 0 then
							AppendTo(strtex, "}");
						fi;
						if j = ll and z = l then
							AppendTo(strtex, " ");
						else
							AppendTo(strtex, " &");
						fi;
					od;
					x := x+ll;
				od;
				AppendTo(strtex, " \\\\ \n");
			od;
			AppendTo(strtex, "\\hline\n\\end{supertabular}\n}\n\n");
			AppendTo(strtex, "\\medskip\n");
		fi;
	fi;
	Print("\nleave compute_fuse()\n\n");
	return Fuse;
end;

sol_factor_fuse := function(Bsol, Bcon, solcon, Fuse, strtex, fv)
# test now which solutions Sol1[1][1][k] (under A itself) factorize 
# over the fusion mapping:
# the result is written into a 3-dim array inv[i][j][k] such that
# the entry is 1, if solution i is invariant under the group Borb[j][k], 
# 0 otherwise
	
	local inv, l, m, n, i, z, ll, j, vl, cand, cand_c, v0, 
	ind, pre, c, pos, len;
	
	Print("\n\nstart sol_factor_fuse():\n");
	if fv > 1 and Length(solcon[1]) > 1 then
		AppendTo(strtex, "Compute the testmatrix $I$ for the fusion process:\n");
		AppendTo(strtex, "\\[ I_{i, j} = \n");
		AppendTo(strtex, "\\left\\{ \\begin{array}{r@{\\quad:\\quad}l}\n");
		AppendTo(strtex, "$k$ & \\text{solution $i$ invariant under group $j$, gives new solution $k$} \\\\ \n");
		AppendTo(strtex, "0 & \\text{else}\n");
		AppendTo(strtex, "\\end{array} \\right. \\]\n\n\\smallskip\n");
	fi;
	
	inv := [];
	l := Length(solcon[1]);
	m := Length(solcon[3][1][1]);
	n := Length(solcon[1][1][1]);
	len := 0;
	for z in [1..l] do
		ind := Bsol[1][z];
		len := len + Length(Bcon[2][ind]);
	od;
	if l > 1 then
		if fv > 1 then
			if len < 10 and len > 4 then
				AppendTo(strtex, "\n\\begin{multicols}{2}\n");
			elif len < 5 then
				AppendTo(strtex, "\n\\begin{multicols}{3}\n");
			fi;
			AppendTo(strtex, "{\\tiny\n\\begin{supertabular}{|");
		fi;
		for z in [1..l] do
			ind := Bsol[1][z];
			if fv > 1 then
				AppendTo(strtex, "*{");
				AppendTo(strtex, String(Length(Bcon[2][ind])));
				AppendTo(strtex, "}{c}");
				AppendTo(strtex, "|");
			fi;
		od;
		if fv > 1 then
			AppendTo(strtex, "}\\hline\n");
		fi;
	fi;
	for i in [1..m] do
		Add(inv, []);
	od;
	for z in [1..l] do
		ll := Length(solcon[1][z]);
		for i in [1..m] do
			inv[i][z]:= [];
		od;
		for j in [1..ll] do
			# number of orbits of Rep1[z][j]
			vl := Length(solcon[1][z][j]);	# number of orbits of Rep1[z][j]
			# Print("group[", z, "][", j, "]:\n");
			# Print("number of orbits on k-subsets: ", vl, "\n");
			cand := [];
			cand_c := [];
			for i in [1..m] do
				Add(cand, 1);
				Add(cand_c, []);
			od;
			for v0 in [1..vl] do
				# compute the pre-image of orbit v0 (canonical ordering)
				pre := [];
				for i in [1..n] do
					if Fuse[z][j][i] = v0 then 	
						# orbit i of A is mapped to v0
						# so pre is an indexset of orbits of A
						Add(pre, i);
					fi;
				od;
				Print("v0=", v0, " pre=", pre, "\n");
				for i in [1..m] do
					if cand[i] = 1 then
						c := is_constant(solcon[3][1][1][i], pre);
						if c = -1 then
							cand[i] := 0;
						else
							Add(cand_c[i], c);
						fi;
					fi;
				od;	# next i
			od; 	# next v0
			# cand_c[i] has maximal length vl (when the solution is constant
			# on the preimages of all orbits under Borb[z][j])
			
			# get the solutions which are invariant under Borb[z][j]:
			for i in [1..m] do
				if cand[i] = 1 then
					pos := find_unsorted(solcon[3][z][j], cand_c[i]);
					if pos = -1 then
						Print("ERROR: solution not found\n");
					else
						# Print("i=", i, " \ncand_c=", cand_c[i], " \nfound at ", pos, "\n");
						Add(inv[i][z], pos);
					fi;
				else
					Add(inv[i][z], 0);
				fi;
			od;	# next i
		od;	# next j
	od;	# next z
	if l > 1 and fv > 1 then
		for i in [1..m] do
			for z in [1..l] do
				ll := Length(solcon[1][z]);
				for j in [1..ll] do
					#AppendTo(strtex, "$");
					if inv[i][z][j] <> 0 then
						AppendTo(strtex, "{\\bf ");
					fi;
					AppendTo(strtex, inv[i][z][j]);
					if inv[i][z][j] <> 0 then
						AppendTo(strtex, "}");
					fi;
					if j = ll and z = l then
						AppendTo(strtex, " ");
					else
						AppendTo(strtex, " &");
					fi;
				od;
			od;
			AppendTo(strtex, " \\\\ \n");
		od;
		AppendTo(strtex, "\\hline\n\\end{supertabular}\n}\n");
		if len < 10 then
			AppendTo(strtex, "\\end{multicols}\n\n");
		fi;
		AppendTo(strtex, "\\medskip\n");
	fi;
	Print("\nleave sol_factor_fuse()\n\n");
	return inv;
end;


sol_true_stab := function(solcon, inv, strtex, fv) 
# get solutions only stabilized by A
	
	local true_stab, i, m, test, z, ll, j, a, S, true_sol, s, n, res, l;
	
	Print("\n\nstart sol_true_stab():\n");
	res := [];
	true_stab := [];
	m := Length(solcon[3][1][1]);
	for i in [1..m] do
	  	Add( true_stab, 1);
	od;
	l := Length(solcon[1]);
	for i in [1..m] do
		test := 0;
		for z in [2..l] do
			ll := Length(solcon[1][z]);
		  	for j in [1..ll] do
		    		a := inv[i][z][j];
				test := test + a;
		  	od; # next j
		od; # next z
		if test <> 0 then
			true_stab[i] := 0;
		fi;
	od; # next i
	
	S := [];
	for i in [1..m] do
	  	if true_stab[i] = 1 then
	    		Add(S, i);
	  	fi;
	od;
	
	if fv > 2 then
		if Length(S)=0 then
			AppendTo(strtex, "All designs have a bigger automorphism group than $A$.\n\n\\smallskip\n");
		elif Length(S)=m then
			AppendTo(strtex, "None of the overgroups arises as automorphism group");
			AppendTo(strtex, "of a design.\n\n\\smallskip\n");
		else		
			AppendTo(strtex, "The solutions, which are not stabilized by one of the above contructed\n");
			AppendTo(strtex, "overgroups of $A$, are the following:\n\n\\smallskip\n");
			AppendTo(strtex, "\\begin{multline*}\nS = (");
			for i in [1..Length(S)] do
				AppendTo(strtex, String(S[i]));
				if i < Length(S) then
					if i mod 25 = 0 then
						AppendTo(strtex, ", \\\\ \n");
					else
						AppendTo(strtex, ", ");
					fi;
				fi;
			od;
			AppendTo(strtex, ")\n\\end{multline*}\n\n\\smallskip\n");
		fi;
	fi;
	
	true_sol := [];
	for i in [1..Length(S)] do
	  	s := solcon[3][1][1][S[i]];
	  	Add(true_sol, s);
	od;
	n := Length(true_sol);
	Print(n, " solution vectors\n");
	
	if fv > 1 then
		if n > 0 and n < m then
			AppendTo(strtex, "We receive ");
			AppendTo(strtex, String(n));
			AppendTo(strtex, " solutions of $A$, which are not invariant under \n");
			AppendTo(strtex, "one of the above constructed overgroups.\n\n\\medskip\n");
		fi;
	fi;
	
	Add(res, true_sol);
	Add(res, S);
	Add(res, n);
	Print("leave sol_true_stab()");
	return res;
end;

