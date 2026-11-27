"""Active circulation sources, distinct from retained anatomical guidance."""


def active_blood_parts(atlas):
    """Reject legacy blood accidentally enabled in a BioDigital replacement.

    Native skin, skeleton and organs remain fitting constraints. Source
    replacement of nerves and lymphatics is a separate, future installation.
    """
    hidden=set(atlas.get('defaultHidden',[]))
    parts=[p for p in atlas['parts']if p['system']in ['arterial','venous']
           and not p.get('suppressed')and not p.get('defaultHidden')and p['id']not in hidden]
    if atlas.get('biodigitalReview'):
        legacy=[p['id']for p in parts if not p['id'].startswith('BDM:')]
        if legacy:
            raise ValueError('Legacy circulation is guidance-only in the BioDigital replacement; active duplicates: '+', '.join(legacy[:8]))
    return parts
