# d4.g
#
# Evi Haberberger
# August/September 1999
#

action_N_NGP_A_on_orbits := function(N_NGP_A, canon, param, KM, strtex, fv) 
# Q gives the operation of N_NGP_A on the orbits and hence the resulting solutions:

	local lg, Q, q, qq, q1, km, km1, i, str, l, j; 

	Print("\n\nstart action_on_orbits():\n");
	lg := Length(GeneratorsOfGroup(N_NGP_A));
	if fv > 0 then
		AppendTo(strtex, "\\subsection{Orbits of the action of $N_{N_G(A)}(P)$ on the solutions}\n\n");
		if fv > 2 then
			AppendTo(strtex, "The action of $N_{N_G(A)}(P)$ on the solutions, which are ");
			AppendTo(strtex, "only invariant under $A$, \nis determined by the ");
			AppendTo(strtex, "following permutations in $Q$:\n\n\\smallskip\n");
			AppendTo(strtex, "%{\\tiny \n");
		fi;
	fi;
	Q := [];
	qq := PermList(canon[2][1]);
	for i in [1..lg] do
	  	km := "KM_file_";
		Append(km, String(KM));
	 	Append(km, "_1_t");
		Append(km, String(param[2]));
		Append(km, "_k");
		Append(km, String(param[3]));
		Append(km, ".txt");
		q := action_on_blocks(km, param[3], GeneratorsOfGroup(N_NGP_A)[i]);
		q1 := qq^-1 * q * qq;
		Add(Q, q1);
		str := String(q1);
		if fv > 2 then
			AppendTo(strtex, "\\begin{align*}\n");
			AppendTo(strtex, str);
			AppendTo(strtex, "\n\\end{align*}\n");
		fi;
	od;
	if fv > 2 then
		AppendTo(strtex, "\n%}\n\n\\medskip\n");
	fi;
	Print("\nleave action_on_orbits()\n\n");
	return Q;
end;

orbit_alg_on_true_sol := function(Q, true_sol, solcon, canon, strtex, fv) 	
# compute canonical representatives of the orbits of Q on true_sol

	local T, R, im, orbit, ol, i, m, l, li, images, a, j, b, d, len, 
	no, k, tlen, sol, S, ind, solu, res;

	Print("\n\nstart orbit_alg_on_true_sol():\n");
	if fv > 2 then
		AppendTo(strtex, "The vector $T$ shows us the canonical representatives of the orbits of this \n");
		AppendTo(strtex, "action on the ");
		AppendTo(strtex, "solutions\n in $S$ (the number in $T$ gives the position in \n");
		AppendTo(strtex, "vector $S$):\n\n\\smallskip\n"); 
	fi;
	
	res := [];
	T :=[];		# transversal of solutions
	R := [];
	im := [];
	orbit := [];
	ol := [];
	m := Length(true_sol[1]);	# number of solutions of A itself	
	d := Length(Q);
	for i in [1..m] do
		orbit[i] := -1;
	od;
	for i in [1..m] do
		#Print("solution ", i, ":\n");
		if orbit[i] = -1 then 	# solution i not yet touched
			R := [i]; 
			l := 1; 
			Add(T, i);
			orbit[i] := Length(T);
			li := 0;
			images := [];
			# find out the orbit of solution i:
			Print("orbit of ", i, ": ");
			while l>0 do 
				a := R[1];
				li := li + 1;
				Add(images, a);
				Print("a = ", a, "\n");
				if l > 1 then 
					for j in [2..l] do
						R[j-1] := R[j];
					od;
				fi;
				l := l-1;
				for j in [1..d] do
					Print("j = ", j, "\n");
					b := Permuted(true_sol[1][a], Q[j]);
					Print("b = ", b, "\n");
					no := -1;
					for k in [1..m] do
						if b = true_sol[1][k] then 
							no := k;
						fi; 
					od;
					if no = -1 then
						Print("\n\nerror: solution not found!\n");
					fi;
					if orbit[no] = -1 then
						# Print(no, " ");
						orbit[no] := Length(T);
						l := l + 1;
						R[l] := no;
					fi;
				od; # next j
			od; # end while-loop
			Add(ol, li);
			Add(im, images);
			# Print("\n");
			j := Length(ol);
			Print("orbit[", j, "]: ", T[j], ", length = ", ol[j], "\n");
		fi;
	od; # next i
	tlen := Length(T);
	# Print("orbit distribution: ", orbit, "\n");
	Print("Length of transversal:", tlen, "\n");
	if fv > 2 then
		AppendTo(strtex, "$T = \\\\\n( ");
		for i in [1..tlen] do
			AppendTo(strtex, String(T[i]));
			if i < tlen then
				if i mod 25 = 0 then
					AppendTo(strtex, ", \\\\ \n");
				else
					AppendTo(strtex, ", ");
				fi;
			fi;
		od; 
		AppendTo(strtex, ")$ \n\n\\smallskip\n");
		if fv > 1 then
			AppendTo(strtex, "The corresponding orbit lengths are $ol = \\\\\n( ");
			for i in [1..tlen] do
				AppendTo(strtex, String(ol[i]));
				if i < tlen then
					if i mod 25 = 0 then
						AppendTo(strtex, ", \\\\ \n");
					else
						AppendTo(strtex, ", ");
					fi;
				fi;
			od;
			AppendTo(strtex, ")$ \n\n\\smallskip\nand the resulting orbit distribution of the solutions is: \\\\ \n$( ");
			for i in [1..m] do
				AppendTo(strtex, String(orbit[i]));
				if i < m then
					if i mod 25 = 0 then
						AppendTo(strtex, ", \\\\ \n");
					else
						AppendTo(strtex, ", ");
					fi;
				fi;
			od;
			AppendTo(strtex, ")$ \n\n");
			# AppendTo(strtex, "\\medskip\n");
		fi;
	fi;
	if fv > 0 then
		AppendTo(strtex, "\n\n\\medskip\nIn the end we get ");
		AppendTo(strtex, String(Length(T)));
		AppendTo(strtex, " orbits of $N_{N_G(A)}(P)$ on the solutions.\n\n");
	fi;
	sol := canon[4][1];		# original solutions under A
	S := true_sol[2];
	if fv > 1 then
		AppendTo(strtex, "The corresponding solutions - i.e. designs - ");
		AppendTo(strtex, "w.r.t the original solution numbers (the index\n");
		AppendTo(strtex, "indicates the corresponding orbit length, that ");
		AppendTo(strtex, "means the number of isomorphic designs):\n\n");
		AppendTo(strtex, "{\\tiny \n");
		len := Length(sol[1]);
		if len < 4 then		
			AppendTo(strtex, "\\begin{multicols}{4}\n");
		elif len > 3 and len < 8 then
			AppendTo(strtex, "\\begin{multicols}{3}\n");
		elif len > 7 and len < 13 then
			AppendTo(strtex, "\\begin{multicols}{2}\n");
		else
			AppendTo(strtex, "\\begin{multicols}{1}\n");
		fi;
		
		for i in [1..Length(T)] do
			ind := S[T[i]];
			solu := sol[ind];
			AppendTo(strtex, "$S_{");
			AppendTo(strtex, String(ind));
			AppendTo(strtex, "}: ");
			AppendTo(strtex, String(solu));
			AppendTo(strtex, "_{");
			AppendTo(strtex, String(ol[i]));
			AppendTo(strtex, "}$ \\\\ \n");
		od;
		AppendTo(strtex, "\\end{multicols} \n}\n");
	fi;
	if fv > 0 then
		AppendTo(strtex, "\n\\medskip\n");
	fi;
	Add(res, T);
	Add(res, ol);
	Add(res, images);
	Print("\nleave orbit_alg_on_true_sol()\n\n");
	return res;
end;
	
file_solutions := function(KM, str, desorb) 
# prints a solution file

	local i, tlen;

	Print("\n\nwrite solution file:\n");
	tlen := Length(desorb[1]);
	PrintTo(str, "");
	for i in [1..tlen] do
		Print("Orbit ", i, ":\n");
		Print("canonical representative: ", desorb[1][i], "\n");
		Print("orbit length: ", desorb[2][i], "\n\n");	
		AppendTo(str, "Orbit ", i, ":\n");
		AppendTo(str, "canonical representative: ", desorb[1][i], "\n");
		AppendTo(str, "orbit length: ", desorb[2][i], "\n\n");	
	od;
end;


span_true_designs := function(T, true_sol, solcon, str)
# span the designs, i.e. the index of the chosen orbits resp. the corresponding
# orbit representative of k-subsets

	local i, j, t, ind, sol, des, design, reps, nborb, r, rep0;
	
	Print("\n\nspan the designs and write into file:\n");
	t := Length(T);	
	PrintTo(str, "# file ", str, "\n# created by dode_proto.g\n\n\n");
	for i in [1..t] do
		ind := T[i];
		sol := true_sol[1][ind];
		des := design_orbits(sol);
		# starting point for the function design_orbits is 0, not 1!	
		# so add 1 to each element of des
		design := des + 1;
		reps := [];
		nborb := Length(design);
		for j in [1..nborb] do
			r := design[j];
			rep0 := solcon[1][1][1][r];
			Add(reps, rep0);
		od;
		Print("design ", i, ": \n", design, "\n");
		Print("corresponding orbit representatives:\n", reps, "\n\n");
		AppendTo(str, "design ", i, ": \n", design, "\n");
		AppendTo(str, "corresponding orbit representatives:\n", reps, "\n\n");
	od;
end;


