# 
# i2.g
#
# Evi Haberberger
# December 1999
#


fuse_with_overgroups := function(W, korb, sol, param, lambda, strtex0, fv)
	local no, l, i, f, str0, rep, repperm, Rep, Rep2, V, v, j,
	fuse, c, k, sol0, Sol0, const, L, S, cmd, solind, S0, S1, len, Repperm, 
	truesol, strtex, res;
	
	# read original KM-file: 
	# solutions and orbit representatives on k-subsets
	
	strtex := [];
	Append(strtex, strtex0);
	Append(strtex, ".tex");
	
	Print("\nstart fuse_with_overgroups():\n\n");
	AppendTo(strtex, "\n\n\\subsection{Fusion of the orbits from the ");
	AppendTo(strtex, "KM-file under action of the constructed overgroups}\n\n");
	
#	Print("orbit reps: ", korb, "\n");
# 	Print("solutions: ", sol, "\n");
	if fv > 0 then
		AppendTo(strtex, "First we get the information from the KM-file:\n\n");
		AppendTo(strtex, "\\begin{itemize}\n");
		AppendTo(strtex, "\\item Original orbit reps of KM: $");
		AppendTo(strtex, String(Length(korb)));
		AppendTo(strtex, "$\n\\item Solutions of KM: $");
		AppendTo(strtex, String(Length(sol)));
		AppendTo(strtex, "$\n\\end{itemize}\n\n");
	fi;

	fuse := [];
	Rep := [];
	Repperm := [];
	l := Length(W);
	for i in [1..l] do
		str0 := "iso_";
		Append(str0, String(i));		
		GeneratorsPermGroup(W[i], str0, param[1]);		
		f := fuse_orbits_by_representatives(korb, str0);
#		Print("New reps:", f[2], "\n");
		Add(fuse, f[2]); 	# vector of image indices
		Add(Rep, f[1]);		# canonical orbit reps under bigger group
		Add(Repperm, ListPerm((1)));	
		# orbits and solutions are already canonical -> identity
		# permutation
	od;
	
	# find the index of the orbit reps in orb[5], and test constance
	truesol := [];
	no := Length(sol);
	for i in [1..no] do
		truesol[i] := 1;
	od;
	
	L := [];
	S := [];
	for i in [1..Length(Rep)] do
		V := [];
		l := Length(Rep[i]);
		for j in [1..l] do
			c := [];
			v := find_indices(fuse[i], j);
			for k in [1..Length(sol)] do
				Add(c, is_constant(sol[k], v));
			od;
			Add(V, c);
		od;
		
		Sol0 := [];
		# now construct the solutions which fuse:
		for k in [1..no] do
			sol0 := [];
			const := 0;
			for j in [1..l] do
				if V[j][k] < 0 then	# no fusion
					const := const + 1;
				else
					Add(sol0, V[j][k]);
				fi;
			od;
			if const = 0 then
				Add(Sol0, sol0);
				truesol[k] := 0;
			fi;
		od;
		if Sol0 <> [] then
			Add(L, i);
			Add(S, Sol0);
		else
			cmd := [];
			Append(cmd, "rm iso_");
			Append(cmd, String(i));
			Print(cmd, "\n");
			Exec(cmd);
		fi;
	od;
	
	solind := [];
	for i in [1..Length(S)] do
		S0 := [];
		for j in [1..Length(S[i])] do
			S1 := [];
			for k in [1..Length(S[i][j])] do
				if S[i][j][k] = 1 then
					Add(S1, k);
				fi;
			od;
			Add(S0, S1);
		od;
		Add(solind, S0);
	od;
	
	# append the info to the latex-file:
	if fv > 0 then
		AppendTo(strtex, "The following groups give solutions i.e. designs:\n\n");
		AppendTo(strtex, "\\begin{itemize}\n");
		for i in [1..Length(L)] do
			AppendTo(strtex, "\\item Group ");
			AppendTo(strtex, String(L[i]));
			AppendTo(strtex, ": \n\\begin{itemize}\n\\item ");
			AppendTo(strtex, String(Length(Rep[i])));
			AppendTo(strtex, " orbits on ");
			AppendTo(strtex, param[3]);
			AppendTo(strtex, "-subsets.\\\\\n");
			if fv > 1 then
				AppendTo(strtex, "{\\tiny \n");
				len := Length(Rep[i][1]);
				if len < 4 then		
					AppendTo(strtex, "\\begin{multicols}{4}\n");
				elif len > 3 and len < 8 then
					AppendTo(strtex, "\\begin{multicols}{3}\n");
				elif len > 7 and len < 17 then
					AppendTo(strtex, "\\begin{multicols}{2}\n");
				else
					AppendTo(strtex, "\\begin{multicols}{1}\n");
				fi;
				for j in [1..Length(Rep[i])] do
					AppendTo(strtex, "$O_{");
					AppendTo(strtex, String(j));
					AppendTo(strtex, "}: \\{");
					for k in [1..Length(Rep[i][j])] do
						AppendTo(strtex, String(Rep[i][j][k]));
						if k < Length(Rep[i][j]) then
							AppendTo(strtex, ", ");
						fi;
					od;
					AppendTo(strtex, "\\}$\\\\\n");
				od;
				AppendTo(strtex, "\\end{multicols}\n}\n");
			fi;
			AppendTo(strtex, "\n\\item ");
			AppendTo(strtex, String(Length(S[i])));
			AppendTo(strtex, " solutions.\\\\\n");
			if fv > 1 then
				AppendTo(strtex, "{\\tiny \n");
				len := Length(solind[i][1]);
				if len < 4 then		
					AppendTo(strtex, "\\begin{multicols}{4}\n");
				elif len > 3 and len < 8 then
					AppendTo(strtex, "\\begin{multicols}{3}\n");
				elif len > 7 and len < 17 then
					AppendTo(strtex, "\\begin{multicols}{2}\n");
				else
					AppendTo(strtex, "\\begin{multicols}{1}\n");
				fi;
				for j in [1..Length(solind[i])] do
					AppendTo(strtex, "$S_{");
					AppendTo(strtex, String(j));
					AppendTo(strtex, "}: \\{");
					for k in [1..Length(solind[i][j])] do
						AppendTo(strtex, String(solind[i][j][k]));
						if k < Length(solind[i][j]) then
							AppendTo(strtex, ", ");
						fi;
					od;
					AppendTo(strtex, "\\}$\\\\\n");
				od;
				AppendTo(strtex, "\\end{multicols}\n}\n");
			fi;
			AppendTo(strtex, "\n\\end{itemize}\n\n");
		od;
		AppendTo(strtex, "\n\\end{itemize}\n\n");
	fi;
	Print("\n\nfusion mapping:\n", fuse, "\n\n");
	
	res := [];
	Add(res, L);
	Add(res, Rep);
	Add(res, S);
	Add(res, Repperm);
	Add(res, fuse);
	Add(res, truesol);
	Print("\nleave fuse_by_overgroups()\n\n");
	return res;
end;


sol_of_conj := function(Bcon, f, param, strtex, fv)
	
	local Rep1, Rep1perm, ll, l, i, t0, z, ind, len, sol, H, rep, lrep,
	Rep1a, Rep1aperm, j, p, rep1, jj, str0, fuse, rep2, rep2perm, Sol1, lr, 
	sol1, k0, solu1, s, s1, truesol, res;
	
	Print("\n\nstart sol_of_conj():\n");
	if fv > 1 and Length(f[1]) > 1 then
		AppendTo(strtex, "\\begin{remark}\n");
		AppendTo(strtex, "Conjugate groups give isomorphic designs, \n");
		AppendTo(strtex, "so the elements of a conjugacy class give similar \n");
		AppendTo(strtex, "solutions as the representative.\n");
		AppendTo(strtex, "\\end{remark}\n\n");
	fi;

	res := [];
	Rep1 := [];
	Rep1perm := [];
	truesol := f[6];
	ll := Length(Bcon[1]);	# nb of conj. classes
	l := Length(f[1]);	# nb of conj. classes with solutions
	for i in [1..ll] do
		Print("conj. class ", i, "\n");
		t0 := 0;
		for z in [1..l] do
			if i = f[1][z] then
				t0 := f[1][z];
				ind := z;
			fi;
		od;
		if t0 > 0 then
			len := Length(Bcon[2][t0]);	# nb of conj. groups in class t0
			Print("orbit length: ", len, "\n");
			sol := f[3][ind];		# sol. of group rep
			H := Bcon[1][t0];		# new (conj.) group
			rep  := f[2][ind];		# k-reps w.r.t. group rep
			Print("representatives: \n", rep, "\n\n\n");
			lrep := Length(rep);
			Rep1a := [];
			Rep1aperm := [];
			for j in [1..len] do
				p := Bcon[3][t0][j];	
				# p perm: group rep -> group
				Print("conjugation permutation: \n", p, "\n\n");
				rep1 := [];
				for jj in [1..lrep] do	# apply to k-reps
					Add(rep1, List(rep[jj], i -> i^p));
				od;
				
				# write file of generators:
				str0 := "iso_";
				Append(str0, String(t0));
				Append(str0, "_");
				Append(str0, String(j));
				GeneratorsPermGroup(Bcon[2][t0][j], str0, param[1]);
		
				fuse := fuse_orbits_by_representatives(rep1, str0);
				rep2 := fuse[1];	# new orbit reps (canonical)
				rep2perm := fuse[2];	# corr. perm for canonical orbit reps
				Add(Rep1a, rep2);
				Add(Rep1aperm, rep2perm);
			od;
			Add(Rep1, Rep1a);
			Add(Rep1perm, Rep1aperm);
		fi;
	od;

	# create the solutions corresponding to those groups
	# w.r.t. canonical orbit representatives
	Sol1 := [];
	lr := Length(Rep1);
	for i in [1..lr] do
	  	sol1 := [];
		for k0 in [1..Length(Rep1perm[i])] do
			Print("k0=", k0, " p=", Rep1perm[i][k0], "\n");
			p := PermList(Rep1perm[i][k0]);
			Print("\nLength Rep1:", lr, "\nLength f[3]:", Length(f[3]), "\n\n");
			Print("solutions: \n", f[3][i], "\n\n");
			sol := f[3][i];
		  	solu1 := [];
	 		for j in [1..Length(sol)] do
	    			s := sol[j];
		    		s1 := Permuted(s, p);
	    			Add(solu1, s1);
			od;
		  	Add(sol1, solu1);
	  	od;
		Add(Sol1, sol1);
	od;
	
	Add(res, Rep1);
	Add(res, Rep1perm);
	Add(res, Sol1);
	Print("\nleave sol_of_conj()\n\n");
	
	return res;
end;

