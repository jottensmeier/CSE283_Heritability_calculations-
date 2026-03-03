# (A) one-time: auth
mkdir -p ~/.ssh && chmod 700 ~/.ssh
ssh-keygen -t rsa -b 4096 -C "jottensmeier@ucsd.edu" -f ~/.ssh/github_rsa
chmod 600 ~/.ssh/github_rsa

cat ~/.ssh/github_rsa.pub   # copy into GitHub settings

cat > ~/.ssh/config <<'EOF'
Host github.com
    HostName github.com
    User git
    IdentityFile ~/.ssh/github_rsa
    IdentitiesOnly yes
EOF
chmod 600 ~/.ssh/config

ssh -T git@github.com

# (B) repo + push
git clone git@github.com:YOURUSER/YOURREPO.git
cd YOURREPO
mkdir -p scripts
cp /path/to/script.py scripts/
git add -A
git commit -m "Add scripts"
git push origin main