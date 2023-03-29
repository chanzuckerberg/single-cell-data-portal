---
name: New analytic events
about: Used to submit new instrumentation requests in Zenhub
title: "Analytics (new): <fill in description>"
labels: analytics, frontend
assignees: ""
---

### Overview

**Tool**: `<Which tool are these events for, i.e. gene expression, explorer>`

- Feature: `<Which feature in that tool are these events for, i.e. filters>`

**Description**: `<Reason for new events request>`
**Expectations**: `<usage of this feature or functionality>`

### New events to add

`<suggested event name>`: `<information that should be stored in the payload (if applicable)>`

- `<Description of what causes the event to fire>`

Checklist

- [ ] Add events to the custom goals section of both the staging and prod Plausible spaces
- [ ] If rdev environment is created for this new feature, once it has been deployed, ping @ainfeld to check analytics
- [ ] Once feature is deployed to Staging ping @ainfeld to check analytics again
- [ ] Verify events are showing up in the Staging space in Plausible correctly
- [ ] Once feature is deployed to Prod ping @ainfeld
- [ ] Confirm Plausible event names in Prod look the same as Staging and are showing up in the Prod space in Plausible correctly

#### People to cc

**PM**: `<tag PM for the tool/feature>`
**DS**: @ainfeld

Privacy policy: https://cellxgene.cziscience.com/privacy
