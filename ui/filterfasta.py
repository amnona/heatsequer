<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>Dialog</class>
 <widget class="QDialog" name="Dialog">
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>400</width>
    <height>300</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>Dialog</string>
  </property>
  <widget class="QDialogButtonBox" name="buttonBox">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>240</y>
     <width>341</width>
     <height>32</height>
    </rect>
   </property>
   <property name="orientation">
    <enum>Qt::Horizontal</enum>
   </property>
   <property name="standardButtons">
    <set>QDialogButtonBox::Cancel|QDialogButtonBox::Ok</set>
   </property>
  </widget>
  <widget class="QLabel" name="label">
   <property name="geometry">
    <rect>
     <x>20</x>
     <y>20</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Fasta:</string>
   </property>
  </widget>
  <widget class="QLineEdit" name="tFileName">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>20</y>
     <width>171</width>
     <height>21</height>
    </rect>
   </property>
  </widget>
  <widget class="QPushButton" name="bBrowse">
   <property name="geometry">
    <rect>
     <x>230</x>
     <y>10</y>
     <width>115</width>
     <height>41</height>
    </rect>
   </property>
   <property name="text">
    <string>Browse</string>
   </property>
  </widget>
  <widget class="QCheckBox" name="cExclude">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>60</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Exclude</string>
   </property>
  </widget>
  <widget class="QLabel" name="lFSExpSamples_3">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>200</y>
     <width>71</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>New name:</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_2">
   <property name="geometry">
    <rect>
     <x>100</x>
     <y>120</y>
     <width>121</width>
     <height>71</height>
    </rect>
   </property>
   <property name="text">
    <string/>
   </property>
   <property name="pixmap">
    <pixmap>arrow_down.png</pixmap>
   </property>
   <property name="scaledContents">
    <bool>false</bool>
   </property>
  </widget>
  <widget class="QLineEdit" name="tNewName">
   <property name="geometry">
    <rect>
     <x>80</x>
     <y>200</y>
     <width>161</width>
     <height>21</height>
    </rect>
   </property>
  </widget>
  <widget class="QCheckBox" name="cOverwrite">
   <property name="geometry">
    <rect>
     <x>270</x>
     <y>200</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>OverWrite</string>
   </property>
  </widget>
  <widget class="QCheckBox" name="cPartialMatch">
   <property name="geometry">
    <rect>
     <x>250</x>
     <y>60</y>
     <width>101</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>PartialMatch</string>
   </property>
  </widget>
 </widget>
 <resources/>
 <connections>
  <connection>
   <sender>buttonBox</sender>
   <signal>accepted()</signal>
   <receiver>Dialog</receiver>
   <slot>accept()</slot>
   <hints>
    <hint type="sourcelabel">
     <x>248</x>
     <y>254</y>
    </hint>
    <hint type="destinationlabel">
     <x>157</x>
     <y>274</y>
    </hint>
   </hints>
  </connection>
  <connection>
   <sender>buttonBox</sender>
   <signal>rejected()</signal>
   <receiver>Dialog</receiver>
   <slot>reject()</slot>
   <hints>
    <hint type="sourcelabel">
     <x>316</x>
     <y>260</y>
    </hint>
    <hint type="destinationlabel">
     <x>286</x>
     <y>274</y>
    </hint>
   </hints>
  </connection>
 </connections>
</ui>
