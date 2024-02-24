Create a virtual environment for the project: 
```shell
python -m venv venv 
```

Start the venv(windows): 
```shell
source venv/Scripts/activate
```

or linux:

```shell
source venv/bin/activate
```

Install the required modules:

```shell
pip install -r requirements.txt
```

The main file is the dna_fountain.py. It can be called with 

```python
python dna_fountain.py
```

Another way is to make it an executable 

```shell
chmod +x dna_fountain.py
```

and add a shebang at the top of the file. That way the file can be executed as following:

```shell
./dna_fountain.py
```