<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>Dialog</class>
 <widget class="QDialog" name="Dialog">
  <property name="windowModality">
   <enum>Qt::WindowModal</enum>
  </property>
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>387</width>
    <height>484</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>Dialog</string>
  </property>
  <widget class="QDialogButtonBox" name="buttonBox">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>450</y>
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
  <widget class="QListWidget" name="lSamples">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>70</y>
     <width>341</width>
     <height>161</height>
    </rect>
   </property>
  </widget>
  <widget class="QListWidget" name="lBacteria">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>260</y>
     <width>341</width>
     <height>181</height>
    </rect>
   </property>
  </widget>
  <widget class="QLabel" name="label">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>10</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Study:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lStudy">
   <property name="geometry">
    <rect>
     <x>80</x>
     <y>10</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_3">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>50</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Samples:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lNumSamples">
   <property name="geometry">
    <rect>
     <x>90</x>
     <y>50</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QLabel" name="lNumBacteria">
   <property name="geometry">
    <rect>
     <x>90</x>
     <y>240</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>?</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_6">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>240</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Bacteria:</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bBiCluster">
   <property name="geometry">
    <rect>
     <x>160</x>
     <y>20</y>
     <width>115</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>BiCluster</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bView">
   <property name="geometry">
    <rect>
     <x>270</x>
     <y>20</y>
     <width>115</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>View</string>
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
