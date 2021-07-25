#!/usr/bin/env bash

# FROM GIST #############################################
#                                                       #
# https://gist.github.com/brenns10/f48e1021e8befd2221a2 #
#                                                       #
#########################################################

# Push HTML files to gh-pages automatically.

# Fill this out with the correct org/repo
ORG=PhyloSofS-Team
REPO=thoraxe
# This probably should match an email for one of your users.
EMAIL=diegozea@gmail.com

set -e

echo "Clone the gh-pages branch outside of the repo and cd into it."
cd ..
git clone -b gh-pages "https://$GH_TOKEN@github.com/$ORG/$REPO.git" gh-pages
cd gh-pages

echo "Update git configuration so I can push."
if [ "$1" != "dry" ]; then
    # Update git config.
    git config user.name "Travis Builder"
    git config user.email "$EMAIL"
fi

echo "Copy in the HTML.  You may want to change this with your documentation path."
cp -R ../$REPO/docs/build/html/* ./

echo "Add and commit changes."
git add -A .
git commit -m "[ci skip] Autodoc commit for $COMMIT."
if [ "$1" != "dry" ]; then
    # -q is very important, otherwise you leak your GH_TOKEN
    git push -q origin gh-pages
fi
