# Happy Remote Development Environment

#### Note: remote dev is still somewhat experimental. Work with #team-happy if you run into trouble.

## Remote Dev Pre-requisites
1. Ensure your `awscli` is configured with the
   [required credentials and profiles](../docs/awscli.md).
1. Make sure you have the *latest version* of the AWS CLI installed. `brew upgrade awscli` if you're not sure:
   ```
   % aws --version
   aws-cli/2.1.8 Python/3.9.0 Darwin/19.6.0 source/x86_64 prompt/off
   ```
1. Run the following command to install prerequisites:
   ```
   pip install -r .happy/requirements.txt
   ```

### Overview
The remote development environment uses a `DEPLOYMENT_STAGE` called `rdev`, but each engineer can run as many remote development *stacks* as they like. Each stack can represent a feature branch, experiment, or whatever's useful to you. Stacks are managed using the remote dev cli utility called `happy`.

The general remote dev workflow is:

1. Make some code changes
1. Run `./scripts/happy create <your-stack-name>` to create a new stack
1. Visit the URL printed by the create step, share it with the team, etc.
1. Run `./scripts/happy logs <your-stack-name> backend` to tail the logs of the data portal api.
1. Make some more code changes
1. Run `./scripts/happy update <your-stack-name>` to update the remote stack with your latest changes.
1. When you don't need your stack anymore, run `./scripts/happy delete <your-stack-name>` to free up remote dev resources.

If you forget which stacks you've created, just run `./scripts/happy list` at any time to list the current remote dev stacks.

### General CLI Usage
The CLI utility is evolving rapidly, so the best reference for which commands are available and how to use them is the CLI itself. All commands support a `--help` flag to print usage docs. For example:

```
% ./scripts/happy create --help
Usage: happy create [OPTIONS] STACK_NAME

  Create a dev stack with a given tag

Options:
  --tag TEXT          Tag name for docker image. Leave empty to generate one
                      automatically.

  --wait / --no-wait  wait for this to complete
  --help              Show this message and exit.
```

### GitHub Action Ingegration
A new stack can also be deployed to remote development environment through GitHub Action integration. Pushing any branch prefixed with "rdev-" will trigger the GH Action workflow to create or update a dev stack, with the stack name equals the part of branch name following the prefix, e.g. pushing branch "rdev-my-dev-branch" will deploy the stack "my-dev-branch" in the remote dev enviroment. This is useful in situation where local connections is slow.

### Authentication
The backend of Happy uses CZI's deployment of Terraform Enterprise (TFE) to deploy and track the resources
for the stacks. This requires logging into TFE to get a long-lived token tied to your user.
The first time you run the happy application, the prompt will give you a command to run to get a token;
follow its prompts.

The long-lived token's access to Terraform Enterprise will periodically expire. The happy CLI will let
you know when this happens, and give you instruction to access the TFE website in your
browser. Loading any TFE web page will reauthorize your token, and you can then re-run your command.

### Warnings

1. Stack name needs to be a valid DNS prefix: starts with a letter, only includes letters, numbers, and dashes, less than 64 characters in length.
1. Yes, you have access to manipulate your teammates' remote dev stacks. This is intentional, to enable collaboration on features. Please use responsibly.

