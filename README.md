# Jupyter book template for Julia Jupyter notebooks

- Canonical: [GitLab](https://gitlab.com/wwtemplates/julia-jupyterbook)
- Mirror: [GitHub](https://github.com/sosiristseng/template-julia-jupyterbook)

## Features

- [Jupyter book](https://jupyterbook.org/index.html) builds `md` and `ipynb` files into a website.
- GitHub actions and GitLab CI/CD build and publish the website whenever changes are committed.
  - The notebook execution results are cached so you can push notebooks with output cell cleared and enjoy the results once the build action is completed.
- Periodically updating Julia dependencies and make a PR if notebooks are executed successfully.
  - For GitHub: you need a pair of SSH keys. (Public key: Deploy key; private key : `SSH_PRIVATE_KEY` actions secret)
  - For GitLab: you need a `GIT_PUSH_TOKEN` [CI/CD variable](https://docs.gitlab.com/ee/ci/variables/index.html), which is a PAT with `write_repo` access.
- Docker images with Julia and jupyter-book dependencies are built automatically by the CI/CD pipelines.
