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
  <widget class="QCheckBox" name="cMitochondria">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>30</y>
     <width>111</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Mitochondria</string>
   </property>
   <property name="checked">
    <bool>true</bool>
   </property>
  </widget>
  <widget class="QCheckBox" name="cChloroplast">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>50</y>
     <width>111</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Chloroplast</string>
   </property>
   <property name="checked">
    <bool>true</bool>
   </property>
  </widget>
  <widget class="QCheckBox" name="cUnknown">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>70</y>
     <width>111</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Unknown</string>
   </property>
   <property name="checked">
    <bool>true</bool>
   </property>
  </widget>
  <widget class="QCheckBox" name="cBacteria">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>90</y>
     <width>141</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Unknown Bacteria</string>
   </property>
   <property name="checked">
    <bool>true</bool>
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
