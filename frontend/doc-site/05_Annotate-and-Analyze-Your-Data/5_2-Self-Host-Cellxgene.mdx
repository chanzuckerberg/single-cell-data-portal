# Self-Host CellxGene


## Hosting cellxgene desktop on the web

Cellxgene desktop is intended to be used by researchers on their local machines. In cases where collaborators are not comfortable installing cellxgene, it is possible to host cellxgene, so your collaborators only need to open the URL you send them.

We don't officially support web deployment, but we've offered some guidance in the following sections on ways to deploy cellxgene to the web.


## General notes and cautions

Please consider the following when deploying cellxgene in any "hosted" environment, especially where access from the broader Internet is possible:



* Information security requires careful configuration of the host environment, including firewall, logging, etc - please follow best practices
* Cellxgene includes features which may be inappropriate for a hosted deployment - you may wish to use the following command line option:` `--disable-diffexp``
*  ``cellxgene launch`` currently uses Flask's development server, which is not recommended for hosted deployment (see [the Flask documentation](https://flask.palletsprojects.com/en/1.1.x/tutorial/deploy/#run-with-a-production-server))
* We have no testing or official support for deployments where multiple users are accessing the same cellxgene instance
* Your cellxgene instance is likely to hang or crash if too many people access it at the same time, especially if they using functions that call the Python backend (such as differential expression, noted above)
* Cellxgene only supports one instance per dataset

If you believe you have found a security-related issue with cellxgene, please report the issue immediately to [[security@chanzuckerberg.com](security@chanzuckerberg.com](mailto:security@chanzuckerberg.com)]


## Configuration Options

The following configuration options require special consideration in any multi-user or hosted environment:

`--disable-diffexp`: the differential expression computation can be resource intensive, in particular for large datasets. If many differential expression calculation requests are made in rapid sequence, it may cause the server CPU or memory resources to be exhausted, and impact the ability of other users to access data. This command line option will disable the differential expression feature, including the removal of the `Differential expression` button.

`--disable-annotations`: annotations, which is enabled by default, may not be appropriate for hosted environments. It will write to the local file system, and in extreme cases could be used to abuse (or exceed) file system capacity on the hosting server. We recommend disabling this with this flag.

`--annotations-file`: this specifies a single file for all end-user annotations, and is incompatible with hosted or multi-user use of cellxgene. Using it will cause loss of user annotation data (ie, the CSV file will be overwritten). If you wish to explore using the annotations feature in a multi-user environment, please refer to the [annotations documentation](LINK TO THE ANNOTATION PAGE ON DOCS SITE), and in particular the `--annotations-dir` flag.


## Community software projects

There are a number of teams building tools or infrastructure to better utilize cellxgene in a multiple user environment. While we do not endorse any particular solution, you may find the following helpful.



* [Novartis Cellxgene Gateway](https://github.com/Novartis/cellxgene-gateway) - a multiple-user and multiple-dataset gateway for cellxgene
* [Interactive Environment in the Galaxy Project](https://galaxyproject.org/) ([patch notes](https://docs.galaxyproject.org/en/release_19.05/releases/19.05_announce.html))

If you know of other solutions, drop us a note and we'll add to this list.


## Live Examples of self-hosted Cellxgene Desktop

Several groups have independently deployed various versions of cellxgene to the web. Check out the cool data that have been hosted using cellxgene!

[Kidney cell atlas](https://www.kidneycellatlas.org/)

[Tabula muris senis](https://tabula-muris-senis.ds.czbiohub.org/)

[Hemocytes](https://hemocytes.cellgeni.sanger.ac.uk/)

[Melanoma](https://melanoma.cellgeni.sanger.ac.uk/)

[Cellxgene Data Portal](https://cellxgene.cziscience.com/)


## Cellxgene Gateway/Apache2 Reverse Proxy with Authentication ▼

Contributors:Fabian Rost and Alexandre Mestiashvili

Hosting Use Case:

- Private self-hosting of multiple datasets for multiple groups

- Groups only have access to datasets that they "own"

- Password Authentication

General description:

Start multiple cellxgene-gateway instances for each of the different user groups and use a reverse proxy for authentication and forwarding to the different cellxgene instances.

Components:

[CellxGene Gateway](https://github.com/Novartis/cellxgene-gateway): the application to be run

[Apache Reverse Proxy](https://httpd.apache.org/docs/2.4/howto/reverse_proxy.html): used as a reverse proxy to redirect web requests and provide authentication

[Tutorial and Code}(https://github.com/mestia/cellxgene-gateway-proxy-example)


## AWS Elastic Bean Stalk ▼

Contributors: The cellxgene team 🧬

Hosting Use Case:



* Cellxgene desktop deployed with elastic beanstalk
* Used to host data in a read-only fashion (can't use annotations)
* Solution for an AWS compute environment

Components:

[CellxGene Desktop](https://github.com/chanzuckerberg/cellxgene): the application to be run

[Elastic Beanstalk](https://aws.amazon.com/elasticbeanstalk/): AWS based solution for deploying and scaling web applications

Note:

This method requires familiarity with

- [IAM](https://aws.amazon.com/iam/)

- [S3](https://aws.amazon.com/s3/)

[Tutorial and Code](([https://github.com/chanzuckerberg/cellxgene/tree/main/backend/czi:hosted/eb](https://github.com/chanzuckerberg/cellxgene/tree/main/backend/czi:hosted/eb))  [THIS IS SUPPOSED TO LINK TO THIS ADDRESS BUT THE LINK IS BROKEN:


## Heroku ▼

Self-hosting with Heroku

You can create a Heroku app via [our provided Dockerfile]([https://github.com/chanzuckerberg/cellxgene/blob/main/Dockerfile](https://github.com/chanzuckerberg/cellxgene/blob/main/Dockerfile)) and [Heroku’s documentation](https://devcenter.heroku.com/articles/build-docker-images-heroku-yml)

You may have to tweak the Docker file like so:

```
1 FROM ubuntu:bionic
2
3 ENV LC_ALL=C.UTF-8
4 ENV LANG=C.UTF-8
5
6 RUN apt-get update && \
7     apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests && \
8     pip3 install cellxgene
9
10 # ENTRYPOINT ["cellxgene"]  # Heroku doesn't work well with ENTRYPOINT
```

and provide a `heroku.yml` file similar to this:

```
1 build:
2   docker:
3     web: Dockerfile
4  run:
5   web:
6    command:
7     - cellxgene launch --host 0.0.0.0 --port $PORT $DATASET # the DATATSET config var must be defined in your dashboard settings.
```

### What is Heroku?

Heroku is a quick and easy way to host applications on the cloud. A Heroku deployment of cellxgene means that the app is not running on your local machine. Instead, the app is installed, configured, and run on the Heroku servers (read: cloud).

On Heroku's servers, applications run on a dyno which are Heroku's implementation and abstraction of containers.

Heroku is one of many options available for hosting instances of cellxgene on the web. Some other options include: Amazon Web Services, Google Cloud Platform, Digital Ocean, and Microsoft Azure.


### Why use Heroku to deploy CellxGene?

What Heroku enables is a quick, non-technical method of setting up a cellxgene instance. No command line knowledge needed. This also allows machines to access the instance via the internet, so sharing a visualized dataset is as simple as sharing a link.

Because cellxgene currently heavily relies on its Python backend for providing the viewer with the necessary data and tooling, it is currently not possible to host cellxgene as a static webpage.

This is a good option if you want to quickly deploy an instance of cellxgene to the web. Heroku deployments are free for small datasets up to around 250MBs in size. See below regarding larger datasets.


### When should I not deploy with Heroku?

The default free dyno offered by Heroku is limited in memory to 512 MBs

- The amount of memory needed for the dyno is roughly the same size as the h5ad file

- Heroku offers tiered paid dynos. More can be found on the [Heroku pricing page​](https://www.heroku.com/pricing)

- Note that this can get _very_ expensive for larger datasets ($25+ a month)

On the free dyno, after 30 minutes of inactivity, Heroku will put your app into a hibernation mode. On the next access, Heroku will need time to boot the dyno back online.

Having multiple simultaneous users requires more memory. This means that the free container size is easily overwhelmed by multiple users, even with small datasets; this can be addressed by purchasing a larger container size

For this facilitated Heroku deployment to work, your dataset must be hosted on a publicly accessible URL

By default, Heroku publically shares your instance to anyone with the URL.



* There are many ways of securing your instance. One quick and simple way is by installing [wwwhisper](https://elements.heroku.com/addons/wwwhisper), a Heroku addon


## AWS/Linux server with basic authentication ▼

Contributors: Lisa Sikkema and Thomas Walzthoeni

Hosting Use Case:

- Private self-hosting on an EC2 instance

- Access by external collaborators and potentially password authentication
