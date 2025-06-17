<ScreeningCriteria version="2025-06-04">
  <Domain name="Population">  
    <Include>  
      <Rule id="POP1">Participants are <emphasis>currently homeless</emphasis> (roofless, houseless, emergency or other temporary accommodation, sofa-surfing).</Rule>  
    </Include>  
    <Exclude>  
      <Rule id="POPX1">Participants are solely <emphasis>formerly</emphasis> homeless.</Rule>  
      <Rule id="POPX2">Participants are residents of asylum-seeker / direct-provision centres.</Rule>  
    </Exclude>  
    <Notes>ETHOS Light categories 1–4 only. If homelessness status is unclear, default to inclusion.</Notes>  
  </Domain>  <!-- ============ DOMAIN 2 : KEY INFORMANTS ============ -->  
  <Domain name="KeyInformants">  
    <Include>  
      <Rule id="KI1">Studies that collect data <emphasis>only</emphasis> from professionals / volunteers reporting on the homeless population.</Rule>  
    </Include>  
    <Notes>No minimum proportion of informants; they must work with homeless people.</Notes>  
  </Domain>  <!-- ============ DOMAIN 3 : GEOGRAPHICAL SETTING ============ -->  
  <Domain name="Setting">  
    <Include>  
      <Rule id="SET1">ata collected in the Republic of Ireland.</Rule>  
      <Rule id="SET2">Multinational studies that include Irish data (even if not broken out separately).</Rule>  
    </Include>  
    <Exclude>  
      <Rule id="SETX1">Records with <emphasis>no</emphasis> Irish data (author affiliation alone is insufficient).</Rule>  
    </Exclude>  
  </Domain>  <!-- ============ DOMAIN 4 : STUDY DESIGN ============ -->  
  <Domain name="StudyDesign">  
    <Include>  
      <Rule id="DES1">Empirical primary or secondary analyses (cross-sectional, cohort, ecological/time-trend, administrative-data reports).</Rule>  
      <Rule id="DES2">Modelling studies that analyse primary or secondary empirical data.</Rule>  
    </Include>  
    <Exclude>  
      <Rule id="DESX1">Single-patient case reports / case series.</Rule>  
      <Rule id="DESX2">Commentaries, editorials, protocols, opinion pieces.</Rule>  
      <Rule id="DESX3">Modelling studies that do <emphasis>not</emphasis> analyse any empirical data.</Rule>  
    </Exclude>  
  </Domain>  <!-- ============ DOMAIN 5 : PUBLICATION TYPE ============ -->  
  <Domain name="PublicationType">  
    <Include>  
      <Rule id="PUB1">Peer-reviewed journal articles (including “online first”).</Rule>  
      <Rule id="PUB2">Conference abstracts.</Rule>  
      <Rule id="PUB3">Systematic reviews <flag code="SR">Include for later mining of primary studies.</flag></Rule>  
    </Include>  
    <Exclude>  
      <Rule id="PUBX1">Pre-prints (e.g. medRxiv).</Rule>  
      <Rule id="PUBX2">Government / NGO reports, white papers, theses.</Rule>  
    </Exclude>  
  </Domain>  <!-- ============ DOMAIN 6 : TOPIC FOCUS ============ -->  
  <Domain name="TopicFocus">  
    <Include>  
      <Rule id="TOP1">Any exposure or intervention with <emphasis>health-related outcomes</emphasis> among the homeless population (e.g. housing policy → mental-health scores).</Rule>  
    </Include>  
    <Exclude>  
      <Rule id="TOPX1">Economic or housing policy studies that do <emphasis>not</emphasis> measure any health outcome.</Rule>  
    </Exclude>  
  </Domain>  <!-- ============ DOMAIN 7 : LANGUAGE ============ -->  
  <Domain name="Language">  
    <Include>  
      <Rule id="LAN1">Abstract available in English, even if full text is another language.</Rule>  
    </Include>  
    <Exclude>  
      <Rule id="LANX1">Abstract not in English, or abstract states full text unavailable in English.</Rule>  
    </Exclude>  
  </Domain>  <!-- ============ DOMAIN 8 : PUBLICATION DATE ============ -->  
  <Domain name="PublicationDate">  
    <Include>  
      <Rule id="DATE1">Earliest publication date (online ahead-of-print or print) is 1 January 2012 or later.</Rule>  
    </Include>  
    <Exclude>  
      <Rule id="DATEX1">Earliest publication date before 2012.</Rule>  
    </Exclude>  
    <Notes>Ignore pre-print posting dates when applying this rule.</Notes>  
  </Domain>  <!-- ============ DOMAIN 9 : GENERAL UNCERTAINTY RULE ============ -->  
  <UncertaintyRule>  
    If reviewers are not highly confident that an abstract fails any inclusion rule, mark it for full-text review.  
  </UncertaintyRule>  <!-- ============ OPTIONAL QUICK-REFERENCE CODES ============ -->  
  <ReferenceCodes>  
    <Code id="INC"   description="Include (meets all rules)"/>  
    <Code id="POP"   description="Exclude — Population"/>  
    <Code id="SET"   description="Exclude — Setting"/>  
    <Code id="DES"   description="Exclude — Design"/>  
    <Code id="PUB"   description="Exclude — Publication type"/>  
    <Code id="DATE"  description="Exclude — Publication date"/>  
    <Code id="LAN"   description="Exclude — Language"/>  
    <Code id="SR"    description="Include — Flag (Systematic Review)"/>  
    <Code id="UNS"   description="Unsure — escalate to full-text review"/>  
  </ReferenceCodes>
</ScreeningCriteria>
