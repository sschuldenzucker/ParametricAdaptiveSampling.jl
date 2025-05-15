# Create docs
makedocs:
    cd docs/; just makedocs

# Manually deploy docs to gh-pages. Requires a clean repo (except for docs/build/).
# You probably only wanna call this from main!
deploydocs: makedocs
    git switch gh-pages
    git rm -rf docs
    # We're deploying from docs/ b/c it makes it much easier to control some meta files (e.g. .gitignore)
    cp -r docs/build docs
    git add .
    git commit -m 'update docs'
    git push
    git switch main
