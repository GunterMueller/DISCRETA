

all: $(DISCRETA_HOME)/obj/$(DISCRETA_ARCH) discreta other main ;

# ruckdeschel

$(DISCRETA_HOME)/obj/$(DISCRETA_ARCH): ;
	cp -r $(DISCRETA_HOME)/obj/empty_dir_structure $(DISCRETA_HOME)/obj/$(DISCRETA_ARCH)

discreta: ; cd DISCRETA; $(MAKE)
ruckdeschel: ; cd RUCKDESCHEL; $(MAKE)
other: ; cd OTHER; $(MAKE)
main: ; cd MAIN; $(MAKE)

clean:
	cd DISCRETA; $(MAKE) clean;
	cd OTHER; $(MAKE) clean;
	cd MAIN; $(MAKE) clean;


