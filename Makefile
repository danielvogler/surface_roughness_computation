# Copyright (c) 2015  Stuart D.C. Walsh and Daniel Vogler,
# All rights reserved.
#
# Contact:
#     Stuart D.C. Walsh  stuart.walsh(at)gmail.com
#     Daniel Vogler  vogler.daniel(at)gmail.com

###########################################################################
#
# Override default makefile settings
#
## Create or edit one of the following files:
##   Make.local.defs
##   Make.local.defs.(username) [e.g. Make.local.defs.stuartwalsh]
## in the build directory to change the default compilation behavior
## `Make.local.defs.username' is provided as a template.
#

-include Make.local.defs
-include Make.local.defs.$(shell whoami)
#
###########################################################################

# export all variables
export

# load all surface classes
classes:
	for classes in $(SURFACE_CLASSES); do \
		$(MAKE) -f MakefileSingle.c SURFACE_CLASS=$$classes; \
		done
