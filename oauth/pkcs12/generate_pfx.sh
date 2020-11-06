if [ ! -f certificate.pfx ]; then
    openssl pkcs12 -export -out certificate.pfx -inkey server.key -in server.crt -password pass:
else
    echo "PFX file already exists, skipping"
fi
