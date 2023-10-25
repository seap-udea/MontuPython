# Version
show:
	@-cat montu/version.py

# GitHub
push:
	git commit -am "New commit"
	git push

#Example: make release RELMODE=release VERSION=0.2.0.2 
release:
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)

# Installation
install:
	pip install -e .

# Cleaning
cleancrap:
	@echo "Cleaning crap..."
	@-find . -name "*~" -delete
	@-find . -name "#*#" -delete
	@-find . -name "#*" -delete
	@-find . -name ".#*" -delete
	@-find . -name ".#*#" -delete
	@-find . -name ".DS_Store" -delete
	@-find . -name "Icon*" -delete
	@-find . -name "*.egg-info*" -type d | xargs rm -fr

cleandist:
	@echo "Cleaning dist..."
	@rm -rf dist/*.*

cleanall:cleancrap cleandist
