# AWS/Linux server with basic authentication

## **Contributors**:

Lisa Sikkema and Thomas Walzthoeni

## Hosting Use Case:

* Private self-hosting on an EC2 instance
* Access by external collaborators and potentially password authentication

## Components:

* [Cellxgene Desktop](https://github.com/chanzuckerberg/cellxgene): the application to be run
* [Nginx](https://www.nginx.com/): used as a reverse proxy to redirect web requests to the WSGI
* [uWSGI](https://uwsgi-docs.readthedocs.io/en/latest/): forwards requests from the Nginx web server to our python flask framework
* [Htpasswd](https://httpd.apache.org/docs/2.4/programs/htpasswd.html): used to allow basic authentication for web requests

## References:

For a general reference, you can [refer](https://hackersandslackers.com/deploy-flask-uwsgi-nginx/) to this article on how to deploy a flask application \(a.k.a. Cellxgene\) to an AWS EC2 instance - or a general linux server using Nginx for receiving web requests, uWSGI to redirect those requests to the python application and htpasswd to add a basic authentication layer over the server.

This general guide can be modified with a few components that would be specific to a cellxgene deployment. First, the flask application which is deployed on your server \(found in the ‘prep your project’ section of the linked article\) should be cellxgene with an appropriate wsgi.py file added as described in the article. Second, the Nginx config file needs to be modified, an example config file may look like the following:

```text
# Default route to login
server {
    listen 443 ssl;
    server_name YOURHOSTNAME;

    # SSL
    ssl_certificate /etc/nginx/conf.d/domain.crt;
    ssl_certificate_key /etc/nginx/conf.d/domain.key;

    location /cellxgene_example {
        proxy_set_header Accept-Encoding "";
        proxy_pass                            
http://localhost:5007/
;
        proxy_set_header Host                 $http_host;
        proxy_set_header X-Real-IP            $remote_addr;
        proxy_set_header X-Forwarded-For      $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto    $scheme;
        # replacements
        sub_filter 'href="static' 'href="cellxgene_example/static';
        sub_filter 'src="static' 'src="cellxgene_example/static';
        sub_filter_types *;
        sub_filter_once off;
        # To add basic authentication to v2 use auth_basic setting.
        auth_basic "Registry realm";
        auth_basic_user_file /etc/nginx/conf.d/nginx.htpasswd;
        client_max_body_size 50M;
    }

}

server {
    listen 80 default_server;

    server_name _;

    return 301 
https://$host$request_uri
;
}
```

With respect to above Nginx config specification, one would start the browser using port 5007 and forward a request from [https://YOURHOSTNAME/cellxgene\_example](https://YOURHOSTNAME/cellxgene_example) to the cxg browser \([http://localhost:5007\](http://localhost:5007%29%29.%20What%20is%20also%20added%20is%20basic_auth%20that%20adds%20a%20login%20request%20with%20a%20username%20and%20password.%20The%20nginx.htpasswd%20file%20is%20generated%20using%20htpasswd.%20If%20you%20want%20to%20add%20additional%20browsers,%20this%20is%20done%20by%20just%20adding%20additional%20location%20stanzas%20and%20running%20the%20cxg%20on%20a%20different%20port.%20The%20additional%20location%20stanzas%20just%20need%20to%20be%20adapted%20accordingly%20%28proxy_pass%20and%20sub_filter%20lines\).

**Note regarding storage and memory on the EC2 instance \(if you implement using AWS\):**

Since datasets will be loaded into the memory of your machine, ram costs could potentially be high, especially if you plan to host multiple large datasets. To avoid higher than necessary costs, you can store data on S3 and only pull the data onto the machine’s memory when the cellxgene explorer is launched. This should be compatible with the approach detailed above \(albeit by tweaking a config variable to point to where the data is located\).

