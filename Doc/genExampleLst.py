#!python3

import yaml

# Load yaml file
with open("exampleLst.yml") as f:
    exampleLst = yaml.safe_load(f)


def single_print(obj, file):
    print('<div class="col">', file=file)
    print('<div class="card h-100 nolink">', file=file)
    print('<a href="{}">'.format(obj['page']), file=file)
    print('<div class="card-img-height fill">', file=file)
    print('<img src="{src}" class="card-img-top" alt="{alt}" height="100%">'.format(**obj['img']), file=f)
    print('</div>', file=file)
    print('<div class="card-body">', file=file)
    print('<h5 class="card-title">{}</h5>'.format(obj['name']), file=file)
    if obj['description'] is not None:
        print('<p class="card-text">{}</p>'.format(obj['description']), file=file)
    print('</div>', file=file)
    print('</a>', file=file)
    print('</div>', file=file)
    print('</div>', file=file)



with open('Examples.md', 'w') as f:
    print('# Example\n', file=f)
    print('@htmlonly\n', file=f)
    print('<div class="row row-cols-1 row-cols-md-3 g-4">', file=f)
    for obj in exampleLst:
          single_print(obj, f)

    print('</div>\n', file=f)
    print('@endhtmlonly', file=f)
