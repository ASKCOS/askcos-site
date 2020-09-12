/* Adapted from https://github.com/64robots/vue-js-panel */

Vue.component('jspanel', {
    template: `
        <div v-if="visible">
            <slot />
        </div>
    `,
    props: {
        visible: {
            type: Boolean,
            default: false
        },
        options: {
            type: Object,
            default: () => ({})
        },
    },
    computed: {
        panelOptions() {
            return Object.assign({onclosed: this.close}, this.options)
        }
    },
    watch: {
        visible(isVisible) {
            if (isVisible) {
                this.createPanel()
            }
        }
    },
    mounted() {
        if (this.visible) {
            this.createPanel()
        }
    },
    methods: {
        async createPanel() {
            await this.$nextTick()
            let options = Object.assign(
                {content: this.$slots.default[0].elm},
                this.panelOptions
            )
            jsPanel.create(options)
        },
        close() {
            this.$emit('close')
            this.$emit('update:visible', false)
        }
    }
})
