# DISCRETA makefile
# Anton Betten 1995
#
# OTHER

all : bin db discreta_solvers largeset mckaynr  \
      lllpacking maxclique solvers tabudes makekm view;

# lllenum spread

bin : ; cd ./BIN ; $(MAKE)
db : ; cd ./DB ; $(MAKE)
discreta_solvers : ; cd ./DISCRETA_SOLVERS; $(MAKE)
largeset : ; cd ./LARGESET ; $(MAKE)
lllenum : ; cd ./LLL_ENUM ; $(MAKE)
mckay : ; cd ./MCKAY ; $(MAKE)
mckaynr : ; cd ./MCKAY_NR ; $(MAKE)
spread : ; cd ./SPREAD ; $(MAKE)
lllpacking : ; cd ./LLL_PACKING ; $(MAKE)
maxclique : ; cd ./MAXCLIQUE ; $(MAKE)
solvers : ; cd ./SOLVERS ; $(MAKE)
tabudes : ; cd ./TABUDES ; $(MAKE)
makekm : ; cd ./MAKE_KM ; $(MAKE)
view : ; cd ./VIEW ; $(MAKE)

clean : ;
	cd ./BIN ; $(MAKE) clean
	cd ./DB ; $(MAKE) clean
	cd ./DISCRETA_SOLVERS ; $(MAKE) clean
	cd ./LARGESET ; $(MAKE) clean
	cd ./LLL_ENUM ; $(MAKE) clean
	cd ./MCKAY ; $(MAKE) clean
	cd ./MCKAY_NR ; $(MAKE) clean
	cd ./SPREAD ; $(MAKE) clean
	cd ./LLL_PACKING ; $(MAKE) clean
	cd ./MAXCLIQUE ; $(MAKE) clean
	cd ./SOLVERS ; $(MAKE) clean
	cd ./TABUDES ; $(MAKE) clean
	cd ./MAKE_KM ; $(MAKE) clean
	cd ./VIEW ; $(MAKE) clean


