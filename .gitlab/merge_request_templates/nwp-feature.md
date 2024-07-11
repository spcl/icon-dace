<!--
ICON

---------------------------------------------------------------
Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
Contact information: icon-model.org

See AUTHORS.TXT for a list of authors
See LICENSES/ for license information
SPDX-License-Identifier: CC0-1.0
---------------------------------------------------------------
-->
## Feature Description
<!--Describe your feature and important technical aspects-->

## Compliance Rules
- One feature per merge request!
- Precise description of changeset and scope.
- Follow the [ICON NWP Workflow](https://gitlab.dkrz.de/icon/icon-nwp/-/wikis/ICON-NWP-Workflow)

## How to Test
- **Compile code with your preferred compiler locally first!**
- Use ICON's build wrappers for testing with different compilers.
- Use a single BuildBot builder for inaccessible compilers, i.e `runBB(breeze_gcc)` or `runBB(levante)`
- Use `runBB` to launch all builders only after resolving all minor issues.

For further information please read [How to use BuildBot](https://gitlab.dkrz.de/icon/wiki/-/wikis/How-to-use-the-new-buildbot#gitlab-merge-requests-for-collective-builds) section in our wiki.

If your changes require the generation of new reference data, please refer to the [Reference data creation](https://gitlab.dkrz.de/icon/wiki/-/wikis/How-to-use-the-new-buildbot#reference-data-creation) section in our wiki.

## Merge Checklist
#### Finalize Feature
- [ ] _Author:_ Cleanup and update branch with icon-nwp:master
- [ ] _Author:_ Test your code according to [How to test](#how-to-test)

#### Scientific Review
- [ ] _Author:_ Select a scientific reviewer and add the label ![inScientificReview](https://img.shields.io/badge/-inScientificReview-orange)
- [ ] _Scientific Reviewer:_ Review
- [ ] _Scientific Reviewer:_ Approve the merge request using Gitlab's approve button.
- [ ] _Scientific Reviewer:_ Remove label ![inScientificReview](https://img.shields.io/badge/-inScientificReview-orange)


#### Gatekeeper Review
- [ ] _Author:_ Mark the merge request as ready by removing `Draft:`
- [ ] _Author:_ Add label ![reviewRequested](https://img.shields.io/badge/-reviewRequested-red) to request a gatekeeper review
- [ ] _Gatekeeper:_ Review (label ![inReview](https://img.shields.io/badge/-inReview-yellow) signals that the Merge Request was accepted for review)
- [ ] _Gatekeeper:_ Prior to merging, please remove any boilerplate from the MR description, retaining only the _Feature Description_ to maintain clean commit messages.

**Note:** Updating the reference data, rebasing and monitoring the buildbot-tests stays in the responsibility of the author.
