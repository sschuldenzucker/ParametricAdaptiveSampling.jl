# Create docs
makedocs:
    cd docs/; just makedocs

# Deploy docs to gh-pages. Requires a clean repo.
deploydocs:
    cd docs/; just deploydocs
