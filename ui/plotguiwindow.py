<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>Dialog</class>
 <widget class="QDialog" name="Dialog">
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>986</width>
    <height>597</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>Dialog</string>
  </property>
  <widget class="QLabel" name="label">
   <property name="geometry">
    <rect>
     <x>8</x>
     <y>0</y>
     <width>51</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Sample:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lSample">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>0</y>
     <width>201</width>
     <height>21</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QLabel" name="lTaxonomy">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>50</y>
     <width>231</width>
     <height>31</height>
    </rect>
   </property>
   <property name="font">
    <font>
     <pointsize>8</pointsize>
    </font>
   </property>
   <property name="text">
    <string>?</string>
   </property>
   <property name="textFormat">
    <enum>Qt::PlainText</enum>
   </property>
   <property name="wordWrap">
    <bool>true</bool>
   </property>
  </widget>
  <widget class="QLabel" name="label_2">
   <property name="geometry">
    <rect>
     <x>8</x>
     <y>50</y>
     <width>31</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Tax:</string>
   </property>
  </widget>
  <widget class="QListWidget" name="lCoolDB">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>150</y>
     <width>256</width>
     <height>192</height>
    </rect>
   </property>
  </widget>
  <widget class="QPushButton" name="bGetSequence">
   <property name="geometry">
    <rect>
     <x>0</x>
     <y>110</y>
     <width>115</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Get Sequence</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_3">
   <property name="geometry">
    <rect>
     <x>8</x>
     <y>70</y>
     <width>21</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>ID:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lID">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>70</y>
     <width>251</width>
     <height>21</height>
    </rect>
   </property>
   <property name="focusPolicy">
    <enum>Qt::ClickFocus</enum>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_4">
   <property name="geometry">
    <rect>
     <x>20</x>
     <y>340</y>
     <width>61</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Selection:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lSelection">
   <property name="geometry">
    <rect>
     <x>90</x>
     <y>340</y>
     <width>131</width>
     <height>21</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bExport">
   <property name="geometry">
    <rect>
     <x>0</x>
     <y>370</y>
     <width>71</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Export</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bView">
   <property name="geometry">
    <rect>
     <x>140</x>
     <y>370</y>
     <width>71</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>View</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bSave">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>370</y>
     <width>71</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Save</string>
   </property>
  </widget>
  <widget class="QComboBox" name="cSampleField">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>20</y>
     <width>104</width>
     <height>26</height>
    </rect>
   </property>
  </widget>
  <widget class="QLabel" name="lSampleFieldVal">
   <property name="geometry">
    <rect>
     <x>120</x>
     <y>20</y>
     <width>201</width>
     <height>21</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QTabWidget" name="FigureTab">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>420</y>
     <width>261</width>
     <height>171</height>
    </rect>
   </property>
   <property name="currentIndex">
    <number>0</number>
   </property>
   <widget class="QWidget" name="tab">
    <attribute name="title">
     <string>None</string>
    </attribute>
    <widget class="QListWidget" name="lStudies">
     <property name="geometry">
      <rect>
       <x>0</x>
       <y>10</y>
       <width>251</width>
       <height>131</height>
      </rect>
     </property>
    </widget>
   </widget>
   <widget class="QWidget" name="ontologies">
    <attribute name="title">
     <string>Ontologies</string>
    </attribute>
    <widget class="QComboBox" name="cOntology">
     <property name="geometry">
      <rect>
       <x>110</x>
       <y>20</y>
       <width>131</width>
       <height>21</height>
      </rect>
     </property>
    </widget>
    <widget class="QLabel" name="label_5">
     <property name="geometry">
      <rect>
       <x>0</x>
       <y>20</y>
       <width>101</width>
       <height>16</height>
      </rect>
     </property>
     <property name="text">
      <string>Ontology Name:</string>
     </property>
    </widget>
   </widget>
   <widget class="QWidget" name="lineplot">
    <attribute name="title">
     <string>LinePlot</string>
    </attribute>
    <widget class="QComboBox" name="cPlotXField">
     <property name="geometry">
      <rect>
       <x>70</x>
       <y>20</y>
       <width>131</width>
       <height>26</height>
      </rect>
     </property>
    </widget>
    <widget class="QLabel" name="label_6">
     <property name="geometry">
      <rect>
       <x>20</x>
       <y>20</y>
       <width>41</width>
       <height>16</height>
      </rect>
     </property>
     <property name="text">
      <string>X field:</string>
     </property>
    </widget>
    <widget class="QCheckBox" name="cPlotXNumeric">
     <property name="geometry">
      <rect>
       <x>20</x>
       <y>50</y>
       <width>89</width>
       <height>20</height>
      </rect>
     </property>
     <property name="text">
      <string>X numeric</string>
     </property>
    </widget>
    <widget class="QCheckBox" name="cPlotNormalizeY">
     <property name="geometry">
      <rect>
       <x>20</x>
       <y>80</y>
       <width>101</width>
       <height>20</height>
      </rect>
     </property>
     <property name="text">
      <string>Normalize Y</string>
     </property>
    </widget>
   </widget>
  </widget>
  <widget class="QPushButton" name="bDBSave">
   <property name="geometry">
    <rect>
     <x>200</x>
     <y>340</y>
     <width>71</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>DBSave</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_7">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>400</y>
     <width>71</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>#samples:</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_8">
   <property name="geometry">
    <rect>
     <x>140</x>
     <y>400</y>
     <width>71</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>#studies:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lNumSamples">
   <property name="geometry">
    <rect>
     <x>80</x>
     <y>400</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QLabel" name="lNumStudies">
   <property name="geometry">
    <rect>
     <x>200</x>
     <y>400</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bExpInfo">
   <property name="geometry">
    <rect>
     <x>110</x>
     <y>110</y>
     <width>91</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>ExpInfo</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bSampleInfo">
   <property name="geometry">
    <rect>
     <x>190</x>
     <y>110</y>
     <width>101</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Sample Info</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_9">
   <property name="geometry">
    <rect>
     <x>8</x>
     <y>90</y>
     <width>41</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Reads:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lReads">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>90</y>
     <width>221</width>
     <height>21</height>
    </rect>
   </property>
   <property name="focusPolicy">
    <enum>Qt::ClickFocus</enum>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bEnrich">
   <property name="geometry">
    <rect>
     <x>200</x>
     <y>370</y>
     <width>71</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Enrich</string>
   </property>
  </widget>
 </widget>
 <resources/>
 <connections/>
</ui>
