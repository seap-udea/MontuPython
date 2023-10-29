# Version
show:
	@-cat montu/version.py

status:
	@-git status

# GitHub
push:cleangit
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
	@-find . -name "__pycache__" -type d | xargs rm -fr
	
cleandist:
	@echo "Cleaning dist..."
	@rm -rf dist/*.*

cleanall:cleancrap cleandist

cleangit:
	@-rm -rf .git/index.lock
	@-rm -rf .git/HEAD.lock
	@-rm -rf .git/refs/remotes/origin/main.lock
	@-rm -rf .git/refs/heads/main.lock

readme:
	python3 -m nbconvert README.ipynb --to markdown
