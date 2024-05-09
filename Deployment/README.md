# Deployment playbooks for MercuryDPM


## Prerequisites

```
  $ pip install ansible>=6.3.0
```


## hosts file

Create a file `hosts` listing your worker hosts (refer to
`hosts.example`)

```
[workers]
worker1.your-domain.org ansible_user=ubuntu
worker2.your-domain.org ansible_user=ubuntu
```

You will also need a key for connecting to your hosts.
```
  $ KEY=/path/to/key.pem
```

If including `localhost` then you should also specify the connection
type, to avoid making an SSH connection.

```
localhost ansible_user=yourname ansible_connection=local
```


## To provision workers

The following command provisions a worker instance by installing
dependencies, creating the source and build directories, cloning the
repository and running `make MercuryBase`.

```
  $ ansible-playbook -i hosts --private-key=$KEY provision.yml
```

By default this will clone the 'trunk' repository at
`https://bitbucket.org/mercurydpm/mercurydpm.git` and checkout the
`HEAD` version. However, you can override this by providing the
repository URL and commit as extra vars:

```
  $ ansible-playbook -i hosts --private-key=$KEY provision.yml \
        -e repository="https://bitbucket.org/you/your_fork.git" \
        -e version="fe5264c"
```

The commit can be specified either by its hash or by the name of a
branch or tag.


## Running a driver

The `runjob.yml` playbook can be used to run a specified driver, for example:

```
  $ ansible-playbook -i hosts --private-key=$KEY runjob.yml \
        -e host=hostname
        -e driver=Tutorials/Tutorial1
        -e rundir=/tmp/demo
```

where you must provide the hostname of the worker, the path to the
driver to run (relative to the `Drivers` directory), and the working
directory (to which results will be written).

This will cause the specified driver to be (re)built and then run as a
background job on the remote server. As with the provision playbook you
can also specify the repository and version to build, although the
default will be to build from the `HEAD` of the trunk.

**TODO**

  1. Be able to specify command line arguments to the driver.
  2. Configure different hosts more easily for `runjob.yml`.
