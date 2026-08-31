# Running Bioinformatics Analyses on AWS

> *For Olivia - congratulations on graduating from the Ram Lab @ UW. May this be one less hurdle between you and your next experiment.*
> *And for Dr. Alison Greenlaw - welcome to the cloud. Trading lab-based servers for AWS is a real shift, but the core logic is the same: get your data somewhere, do the work, save it somewhere safe. You've got this.*

If you have sequencing data that is too big for your laptop - maybe a few hundred GB of FASTQ files - and you need somewhere with more horsepower to align reads, run featureCounts, or process an nf-core pipeline. This tutorial walks you through renting that horsepower from AWS, running your analysis, and shutting everything down cleanly (so you don't get a surprise bill).

No prior AWS experience needed. By the end, you should understand not just what to type, but why.

---

## 1. The Big Picture

Think of AWS like renting a fully equipped lab bench by the hour, instead of buying one outright.

- EC2 is the bench itself - a computer you can configure with however much CPU and RAM your analysis needs.
- S3 is the building's archive room - cheap, durable, long-term storage for your raw data and final results. It's not where you compute; it's where you stash things before and after.
- You rent the bench (EC2), pull your samples out of the archive (S3), do the work, put your results back in the archive, then give the bench back so you stop paying for it.

That last step matters more than any other line in this tutorial. An EC2 instance left running over a weekend is the single most common way people accidentally rack up a large AWS bill.

### The workflow, visually

```
 Your Laptop
     │  ssh
     ▼
 EC2 Instance  ◄──────────────┐
     │                        │
     │ 1. pull FASTQs         │ 4. push results
     ▼                        │
   S3 Bucket (raw data)   S3 Bucket (results)
                              ▲
     2. run alignment ────────┘
     3. write results
```

Five things happen in order: fetch data → compute → save results → tear down the instance. Everything below is just the detail behind those four verbs.

---

## 2. Vocabulary You Will Actually Use

Do not memorize this - skim it once, then refer back when a term comes up later.

| Term | What it actually is | Analogy |
|---|---|---|
| EC2 | A virtual computer you rent by the hour/second | The lab bench |
| EBS | The hard drive attached to that computer | The drawer under the bench where you keep working files |
| Snapshot | A backup copy of an EBS drive | A photocopy of everything in that drawer, stored safely elsewhere |
| S3 | Long-term, cheap storage, separate from any running computer | The archive room down the hall |
| Key Pair | A cryptographic key that proves it's you logging in | Your personal key to the building - lose it, and you can't get back in |
| Security Group | A firewall - controls what traffic can reach your instance | The bouncer deciding who's allowed through the door |

---

## 3. Before You Launch: Four Decisions

When you create an EC2 instance, AWS asks you to choose four things. Here's how to think about each one.

### a) Instance type - how much CPU/RAM?

This is the size of your rented bench, and it directly affects your hourly cost. A good way to think about it:

```
T2 Micro        →  for testing the workflow / getting familiar with the console
C3.4XLarge      →  for an actual compute job (alignment, featureCounts, etc.)
```

Start small with a T2 Micro just to confirm everything connects and works the way you expect. Once you're ready to run a real job on real data, switch to something like a C3.4XLarge for the extra CPU and memory. If your pipeline still runs out of memory, the fix is to size up further next time, not to fight with the software.

Or give your file size related prompt to Claude and budget and ask for suitable instance type.

### b) Image (AMI) - which operating system?

This is the software your bench comes pre-loaded with. Amazon Linux is a great default starting point - it comes with a lot of common dependencies pre-installed, which saves you setup time, and it's well-supported within the AWS ecosystem itself. Ubuntu is the other common choice in bioinformatics (most external tools and tutorials assume it), but Amazon Linux is the easier on-ramp if you're just getting started:

```
Amazon Linux  (or Ubuntu 22.04, if your tools specifically expect it)
```

One thing to note if you choose Amazon Linux: your SSH username and package manager will be different from Ubuntu - more on that in Step 4.

### c) Storage (EBS) - how much disk space?

This is your working drawer space for FASTQs, BAMs, intermediate files, and results. Sequencing data is bulky, so don't be stingy:

```
2.4 TB, GP3 volume type
```

GP3 is the default volume type, and it's a good fit here because it's cost-effective and performs well for general-purpose workloads that don't have heavy, constant IO demands - which covers most bioinformatics pipelines. If 2.4 TB feels like overkill for your dataset, scale it down - but as a rule of thumb, budget for roughly 3-4x the size of your raw FASTQs, since alignment produces large intermediate files (sorted BAMs, indices, temp files) before you get to your final results.

### d) Key Pair & Security Group - how do you get in, and who else can?

- Key Pair: generate one and save the `.pem` file somewhere safe. There is no "forgot password" recovery for this - if you lose it, you lose access to that instance.
- Security Group: at minimum, allow inbound traffic on:

```
SSH - Port 22
```

This is the only door you need open for this workflow.

---

## 4. Connecting to Your Instance

Once your instance is running, AWS gives you a public IP address. From your laptop's terminal:

```bash
# Go to wherever you saved the key
cd ~/Downloads

# Lock down the key's permissions (AWS requires this - SSH will refuse to use an unprotected key)
chmod 400 my-key.pem

# Connect - username depends on which AMI you picked
ssh -i my-key.pem ec2-user@<PUBLIC-IP>     # Amazon Linux
ssh -i my-key.pem ubuntu@<PUBLIC-IP>       # Ubuntu
```

You will know it worked when your prompt changes to something like:

```
[ec2-user@ip-172-31-xx-xx ~]$     # Amazon Linux
ubuntu@ip-172-31-xx-xx:~$         # Ubuntu
```

You are now inside the cloud computer, not your laptop. Everything you type from here runs on AWS.

Quick note on package managers: if you went with Amazon Linux, you'll use `yum` (or `dnf` on newer versions) to install software. On Ubuntu, that's `apt`. Same idea, different command.

---

## 5. Organize Before You Work

A few seconds of setup now saves a lot of confusion later, especially once you have multiple samples and log files flying around.

```bash
mkdir aws-bioinformatics-project
cd aws-bioinformatics-project
mkdir data scripts results logs
```

```
aws-bioinformatics-project/
├── data/       ← raw FASTQs go here
├── scripts/    ← your analysis scripts
├── results/    ← anything you intend to keep
└── logs/       ← stdout/stderr from your runs
```

Keeping `results/` and `logs/` separate from `data/` makes Step 8 (uploading) much less error-prone - you just sync two folders, not a mess of mixed files.

---

## 6. Pulling Data Down from S3

First, confirm you can see your bucket:

```bash
aws s3 ls
```

Then copy your sequencing files into the `data/` folder you just made:

```bash
aws s3 sync s3://my-bucket/fastq/ data/
```

`sync` is preferable to `cp` here - if the transfer is interrupted and you run it again, it only copies what's missing instead of starting over.

Double-check the files actually arrived before you go further:

```bash
ls data
```

---

## 7. Running Your Analysis

This is the part that's specific to your actual experiment - common tools include:

- Rsubread / featureCounts - read counting for RNA-seq
- Nextflow + nf-core pipelines - standardized, reproducible pipelines for alignment, variant calling, etc.

While something is running, keep an eye on whether your instance can actually handle it:

```bash
htop      # live CPU and memory usage
free -h   # how much RAM is free, human-readable
df -h     # how much disk space is left, human-readable
```

If `df -h` starts creeping toward 100%, stop and clean up intermediate files before your job crashes mid-run from a full disk - that's a much worse time to discover the problem.

---

## 8. Don't Let a Disconnect Kill Your Job

If your laptop's connection drops - wifi hiccup, your laptop goes to sleep, you close the lid - a job running directly in that SSH session dies with it. `screen` solves this by running your job in a session that keeps going on the AWS side, independent of your laptop's connection.

```bash
# Start a named session
screen -S alignment

# ...run your analysis here as normal...

# Detach without stopping anything (Ctrl+A, then D)
```

Your job keeps running. You can close your laptop, walk away, come back tomorrow, and:

```bash
screen -r alignment
```

...and you are looking at the same session, picking up exactly where you left off.

---

## 9. Saving Your Results Back to S3

Once your analysis finishes, push everything worth keeping back to the archive room - your EC2 instance is temporary, but S3 is durable.

```bash
aws s3 sync results/ s3://my-results/project/results/
aws s3 sync logs/ s3://my-results/project/logs/
```

---

## 10. Shutting Down Safely

This is the step that actually saves you money, so don't skip the verification.

Before stopping or terminating, confirm:
- ✅ Results are uploaded to S3
- ✅ Scripts are backed up somewhere (S3, GitHub, etc.)
- ✅ Logs are saved

Verify it landed, don't just trust that the sync command "probably worked":

```bash
aws s3 ls s3://my-results/project/results/
```

Then choose:

| Action | When to use it |
|---|---|
| Stop | You'll likely come back to this same instance later. EBS storage charges continue, but you stop paying for compute. |
| Terminate | The project is fully done. Everything on the instance (other than what you saved to S3) is gone for good. |

---

## 11. Mistakes That Catch Almost Everyone Once

- ❌ Terminating an instance before confirming results actually made it to S3
- ❌ Losing the `.pem` key (write down where you saved it - today)
- ❌ Running out of disk mid-job because intermediate files weren't accounted for
- ❌ Skipping `screen` and losing hours of progress to a dropped connection
- ❌ Downloading data into the wrong folder, then losing track of what's where
- ❌ Leaving an instance running overnight "just in case" - this is the #1 source of unexpectedly large bills

---

## 12. Quick Reference / Cheat Sheet

```bash
# Connect to your instance
ssh -i my-key.pem ec2-user@<PUBLIC-IP>   # Amazon Linux
ssh -i my-key.pem ubuntu@<PUBLIC-IP>     # Ubuntu

# Pull data down from S3
aws s3 sync s3://bucket/data data/

# Monitor while running
htop
free -h
df -h

# Push results up to S3
aws s3 sync results/ s3://bucket/results/

# Reconnect to a screen session after disconnecting
screen -r
```

---

## One-Sentence Summary

Rent the bench (EC2) → pull samples from the archive (S3) → compute → push results back to the archive (S3) → return the bench (stop/terminate) - and always verify step four before step five.
