# Contributing to Single Cell Data Portal

## Pull Request Protocols

The followings guidelines allows the core developers to manage the workload of performing PR reviews with a 24 hour SLO more manageable.

- Please do not request code reviews for PRs that are failing unit tests, deployment, or lint, unless it was failing prior to your PR (in which case I’d recommend figuring out why it was failing prior to your PR to make sure it is being addressed).

- Please do not request code reviews on draft PRs unless you're specifically asking for someone to scan something or validate (in which case, please contact them beforehand).

- Once you've fully addressed all the comments, please re-request a review and/or offline ping. Sometimes it is difficult figuring out which PRs have fully addressed all the comments and which ones have not.

- Please resolve comments that have been addressed.

- If someone has not responded on a PR for a while, please ping them directly. It's preferable that PRs get closed instead of waiting indefinitely for someone to review while new PRs keep piling up.

- For not urgent but safe/trivial PRs, please use label `bot/merge` and request reviews, so `czimergebot` will automatically merge the PR for you once you get at least one approval + when all Github PR checks pass

- For urgent but safe/trivial PRs, please feel free to force merge now, but still tag reviewers, in case they spot low priority fixes that can be done as a followup PR

- Please title your PR's with the intended squash commit message. This title should follow [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) standard (`<type>[optional scope]: <description>`). PR titles will be linted via the [lint-pr workflow](https://github.com/chanzuckerberg/single-cell-data-portal/blob/main/.github/workflows/lint-pr.yml).
