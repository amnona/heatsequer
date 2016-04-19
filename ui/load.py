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
  <widget class="QPushButton" name="bLoadBrowseTable">
   <property name="geometry">
    <rect>
     <x>240</x>
     <y>20</y>
     <width>115</width>
     <height>41</height>
    </rect>
   </property>
   <property name="text">
    <string>Browse</string>
   </property>
  </widget>
  <widget class="QLineEdit" name="tLoadTable">
   <property name="geometry">
    <rect>
     <x>70</x>
     <y>30</y>
     <width>171</width>
     <height>21</height>
    </rect>
   </property>
  </widget>
  <widget class="QLabel" name="label">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>30</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Table:</string>
   </property>
  </widget>
  <widget class="QLabel" name="label_2">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>70</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Mapping:</string>
   </property>
  </widget>
  <widget class="QLineEdit" name="tLoadMap">
   <property name="geometry">
    <rect>
     <x>70</x>
     <y>70</y>
     <width>171</width>
     <height>21</height>
    </rect>
   </property>
  </widget>
  <widget class="QPushButton" name="bLoadBrowseMap">
   <property name="geometry">
    <rect>
     <x>240</x>
     <y>60</y>
     <width>115</width>
     <height>41</height>
    </rect>
   </property>
   <property name="text">
    <string>Browse</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bLoadLoad">
   <property name="geometry">
    <rect>
     <x>30</x>
     <y>220</y>
     <width>131</width>
     <height>71</height>
    </rect>
   </property>
   <property name="text">
    <string>Load</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bLoadCancel">
   <property name="geometry">
    <rect>
     <x>220</x>
     <y>220</y>
     <width>131</width>
     <height>71</height>
    </rect>
   </property>
   <property name="text">
    <string>Cancel</string>
   </property>
  </widget>
  <widget class="QCheckBox" name="cLoadDeblurred">
   <property name="geometry">
    <rect>
     <x>80</x>
     <y>110</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Deblurred</string>
   </property>
   <property name="iconSize">
    <size>
     <width>20</width>
     <height>20</height>
    </size>
   </property>
   <property name="checked">
    <bool>true</bool>
   </property>
  </widget>
  <widget class="QLineEdit" name="tLoadName">
   <property name="geometry">
    <rect>
     <x>140</x>
     <y>160</y>
     <width>211</width>
     <height>21</height>
    </rect>
   </property>
  </widget>
  <widget class="QLabel" name="label_3">
   <property name="geometry">
    <rect>
     <x>20</x>
     <y>160</y>
     <width>121</width>
     <height>16</height>
    </rect>
   </property>
   <property name="palette">
    <palette>
     <active>
      <colorrole role="WindowText">
       <brush brushstyle="SolidPattern">
        <color alpha="255">
         <red>56</red>
         <green>87</green>
         <blue>231</blue>
        </color>
       </brush>
      </colorrole>
     </active>
     <inactive>
      <colorrole role="WindowText">
       <brush brushstyle="SolidPattern">
        <color alpha="255">
         <red>56</red>
         <green>87</green>
         <blue>231</blue>
        </color>
       </brush>
      </colorrole>
     </inactive>
     <disabled>
      <colorrole role="WindowText">
       <brush brushstyle="SolidPattern">
        <color alpha="255">
         <red>127</red>
         <green>127</green>
         <blue>127</blue>
        </color>
       </brush>
      </colorrole>
     </disabled>
    </palette>
   </property>
   <property name="text">
    <string>Experiment Name:</string>
   </property>
  </widget>
  <widget class="QCheckBox" name="cMetabolite">
   <property name="geometry">
    <rect>
     <x>170</x>
     <y>110</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Metabolite</string>
   </property>
   <property name="iconSize">
    <size>
     <width>20</width>
     <height>20</height>
    </size>
   </property>
   <property name="checked">
    <bool>false</bool>
   </property>
  </widget>
  <widget class="QCheckBox" name="cEMP">
   <property name="geometry">
    <rect>
     <x>270</x>
     <y>110</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="toolTip">
    <string>change map sampleid to lowercase</string>
   </property>
   <property name="text">
    <string>EMP</string>
   </property>
   <property name="iconSize">
    <size>
     <width>20</width>
     <height>20</height>
    </size>
   </property>
   <property name="checked">
    <bool>false</bool>
   </property>
  </widget>
 </widget>
 <resources/>
 <connections/>
</ui>
