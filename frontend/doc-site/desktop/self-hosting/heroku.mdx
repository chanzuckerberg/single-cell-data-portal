# Heroku

## Self-hosting with Heroku

You can create a Heroku app via [our provided Dockerfile here](https://github.com/chanzuckerberg/cellxgene/blob/main/Dockerfile) and [Heroku's documentation](https://devcenter.heroku.com/articles/build-docker-images-heroku-yml).

You may have to tweak the `Dockerfile` like so:

```text
FROM ubuntu:bionic

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && \
    apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests && \
    pip3 install cellxgene

# ENTRYPOINT ["cellxgene"]  # Heroku doesn't work well with ENTRYPOINT
```

and provide a `heroku.yml` file similar to this:

```text
build:
  docker:
    web: Dockerfile
run:
  web:
    command:
      - cellxgene launch --host 0.0.0.0 --port $PORT $DATASET # the DATATSET config var must be defined in your dashboard settings.
```

## What is Heroku?

Heroku is a quick and easy way to host applications on the cloud.

A Heroku deployment of cellxgene means that the app is not running on your local machine. Instead, the app is installed, configured, and run on the Heroku servers \(read: cloud\).

On Heroku's servers, applications run on a [dyno](https://www.heroku.com/dynos) which are Heroku's implementation and abstraction of containers.

Heroku is one of many options available for hosting instances of cellxgene on the web. Some other options include: Amazon Web Services, Google Cloud Platform, Digital Ocean, and Microsoft Azure.

## Why use Heroku to deploy cellxgene?

What Heroku enables is a quick, non-technical method of setting up a cellxgene instance. No command line knowledge needed. This also allows machines to access the instance via the internet, so sharing a visualized dataset is as simple as sharing a link.

Because cellxgene currently heavily relies on its Python backend for providing the viewer with the necessary data and tooling, it is currently not possible to host cellxgene as a static webpage.

This is a good option if you want to quickly deploy an instance of cellxgene to the web. Heroku deployments are free for small datasets up to around 250MBs in size. See below regarding larger datasets.

## When should I not deploy with Heroku?

* The default free dyno offered by Heroku is limited in memory to 512 MBs
  * The amount of memory needed for the dyno is roughly the same size as the h5ad file
  * Heroku offers tiered paid dynos. More can be found on the [Heroku pricing page](https://www.heroku.com/pricing)
  * Note that this can get _very_ expensive for larger datasets \($25+ a month\)
* On the free dyno, after 30 minutes of inactivity, Heroku will put your app into a hibernation mode. On the next access, Heroku will need time to boot the dyno back online.
* Having multiple simultaneous users requires more memory. This means that the free container size is easily overwhelmed by multiple users, even with small datasets; this can be addressed by purchasing a larger container size
* For this facilitated Heroku deployment to work, your dataset must be hosted on a publicly accessible URL
* By default, Heroku publically shares your instance to anyone with the URL.
  * There are many ways of securing your instance. One quick and simple way is by installing [wwwhisper](https://elements.heroku.com/addons/wwwhisper), a Heroku addon

