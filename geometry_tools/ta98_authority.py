"""Explicit corrections to the cached TA98 extract, with primary provenance.

The SQLite cache remains untouched. Both viewer generation and construction
inventories must use the same reviewed parents and sex-qualified identities.
"""
FEMALE_URETHRA = 'https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/09.2.03%20TA98%20EN.htm'
MALE_URETHRA = 'https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/09.4.02%20TA98%20EN.htm'
BLADDER = 'https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/08.3.01%20TA98%20EN.htm'
OPENING_STATES = {
    'A09.2.03.001': ('A08.3.01.028F', 'female', FEMALE_URETHRA),
    'A09.2.03.002': ('A08.3.01.028F', 'female', FEMALE_URETHRA),
    'A09.4.02.001': ('A08.3.01.027M', 'male', MALE_URETHRA),
    'A09.4.02.002': ('A08.3.01.027M', 'male', MALE_URETHRA),
}
# Listed explicitly on the bladder section page, under its canonical pelvic
# muscle parent elsewhere. This is a context occurrence, never a second parent.
CONTEXT_REFERENCES = {'A08.3.01.001': {'A04.5.03.017': BLADDER}}


def corrected_terms(source):
    result={code:dict(row)for code,row in source.items()}
    for code,(parent,sex,url)in OPENING_STATES.items():
        if code not in result:continue
        row=result[code]
        if parent not in result:
            raise ValueError('The sex-qualified opening authority is incomplete: '+parent)
        old=row['parent_id']
        if old not in {'A08.3.01.027M','A08.3.01.028F'}:
            raise ValueError('Unexpected cached opening-state parent requires fresh review: '+code)
        if old!=parent:
            row['source_parent_id']=old
            row['parent_id']=parent
            row['parent_correction']='Sex-qualified opening ancestry verified against the primary TA98 section'
            row['parent_correction_source']=url
            if 'parent_name'in row:row['parent_name']=result[parent]['name_en']
        for field,value in [('female_gender',int(sex=='female')),('male_gender',int(sex=='male'))]:
            if row.get(field)!=value:
                row['source_'+field]=row.get(field)
                row[field]=value
        row['sex_correction_source']=url
    return result


def term_sex(terms,code):
    """Use explicit authority gender or a sex-qualified ancestor, not names."""
    seen=set()
    while code in terms and code not in seen:
        seen.add(code);row=terms[code]
        female=bool(row.get('female_gender'));male=bool(row.get('male_gender'))
        if female!=male:return 'female'if female else'male'
        if code.endswith('F'):return 'female'
        if code.endswith('M'):return 'male'
        code=row.get('parent_id')
    return None
