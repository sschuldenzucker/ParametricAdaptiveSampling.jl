# Create docs
makedocs:
    cd docs/; just makedocs

# Manually deploy docs to gh-pages. Requires a clean repo (except for docs/build/).
# You probably only wanna call this from main!
deploydocs: makedocs
    git switch gh-pages
    git rm -rf .
    cp -r docs/build/* ./
    git add .
    git commit -m 'update docs'
    git push
    git switch main
